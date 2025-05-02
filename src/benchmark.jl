#! /usr/bin/env julia
"""
Run jet reconstruction over a standard set of input files, with
varying initial particle densities
"""

using ArgParse
using Logging
using JSON
using CSV
using DataFrames
using EnumX
using CodecZlib

using LorentzVectorHEP
using JetReconstruction

# Backends for the jet reconstruction
@enumx T=Code Backends JetReconstruction Fastjet AkTPython AkTNumPy
const AllCodeBackends = [String(Symbol(x)) for x in instances(Backends.Code)]

# Parsing for Enum types
function ArgParse.parse_item(opt::Type{E}, s::AbstractString) where {E <: Enum}
    insts = instances(E)
    p = findfirst(x -> Symbol(x) == Symbol(s), insts)
    if isnothing(p)
        throw(ErrorException("Invalid value for enum $opt: $s"))
    end
    return insts[p]
end

function determine_backend_version(code::Backends.Code)
    if code == Backends.JetReconstruction
        return string(VERSION)
    elseif code == Backends.Fastjet
        return "unknown"
    elseif code in (Backends.AkTPython, Backends.AkTNumPy)
        output = read(`python --version`, String)
        m = match(r"Python (\d+\.\d+\.\d+)", output)
        if isnothing(m)
            return "unknown"
        end
        return m[1]
    end
    "unknown"
end

function determine_backend(code::Backends.Code)
    if code == Backends.JetReconstruction
        return "Julia"
    elseif code == Backends.Fastjet
        return "C++"
    elseif code in (Backends.AkTPython, Backends.AkTNumPy)
        return "Python"
    end
    "unknown"
end

function hepmc3gunzip(input_file::AbstractString)
    unpacked_file = replace(input_file, ".gz" => "")
    if !isfile(unpacked_file)
        @info "Unpacking $(input_file) to $(unpacked_file)"
        in = GzipDecompressorStream(open(input_file))
        out = open(unpacked_file, "w")
        write(out, in)
        close(in)
        close(out)
    end
    unpacked_file
end

function ipmi_command()
    output = read(`sudo ipmitool dcmi power reading`, String)
    t = Int64(time_ns())
    s = tryparse(Float64, match(r"Instantaneous power reading:\s+([\d\.]+) Watts", output)[1])
    (s, t)
end

#returns power in uj at current sample time, detecting "loops" since the last measurement by incrementing the loops variable
#you must sample ratp sufficiently frequently that only 1 loop happens per sample period [probably > 0.05Hz]
#(power, loops)  - you will need to calculate total power by multiplying loops by max_ratp_read() and adding to the power
function ratp_read(prev_v, loops)
	new_v = 0.0
	open("/sys/class/powercap/intel-rapl:0/energy_uj", "r") do f
		new_v = tryparse(Float64,read(f, String))
		if new_v < prev_v
			loops += 1
		end
	end
	(new_v, loops)
end

#maximum power range (the value at which ratp_read loops) in uj
#total power from a ratp_read sampled sufficiently often is power*max_ratp_read*loops [see ratp_read]
function max_ratp_read()
	f = open("/sys/class/powercap/intel-rapl:0/max_energy_range_uj", "r") 
	parse(Float64,read(f, String))
end

function ipmi_chan(nsamples::Integer)
    Channel{Float64}(spawn=true) do ch 
        ipmi_time = [(0.0,0.0), (0.0,0.0)] 
	ratp_sample = [0.0,0.0]
        max_ratp = max_ratp_read()
        idx = 1
        energy_uj = 0.0

        rdy() = (isready(ch))
    
        for irun in 1:nsamples
            take!(ch) #block until we tell it to start
            ipmi_time[idx] = ipmi_command()
	    sample_energy, loops = ratp_read(0, 0)
	    ratp_sample[idx] = sample_energy
            idx = 3 - idx #toggle 1 and 2
            while !rdy()
                    timedwait(rdy, 5) #sample every 5 seconds or until we're told to stop
                    ipmi_time[idx] = ipmi_command()
 		    ratp_sample[idx],loops = ratp_read(ratp_sample[3-idx],loops)  #need to compare to previous to discover discontinuities when it overflows
                    idx = 3 - idx
                    #trapezium rule
                    energy_uj += (ipmi_time[1][1]+ipmi_time[2][1])*abs(ipmi_time[2][2]-ipmi_time[1][2])/2.0
            end
            take!(ch) #take the value passed to us to stop our counter 
		#the function we pass to uses us as its timescale, so if we pass uj then uj/us = W, which is what we want
            put!(ch, energy_uj * 1.e-3 ) #uJoules from nj
            put!(ch, (ratp_sample[3-idx] + max_ratp*loops - sample_energy) ) #uj from uj 
            energy_uj = 0.0
        end
    end    
end
    
function julia_jet_process_avg_time(events::Vector{Vector{T}};
    ptmin::Float64 = 5.0,
    radius::Float64 = 0.4,
    p::Union{Real, Nothing} = nothing,
    algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
    strategy::RecoStrategy.Strategy,
    nsamples::Integer = 1,
    repeats::Int = 1, 
    ipmi::Bool = false ) where {T <: JetReconstruction.FourMomentum}
    @info "Will process $(size(events)[1]) events, repeating $(repeats) time(s)"

    # Set consistent algorithm and power
    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
    algorithm = algorithm)
    
    n_events = length(events)
    ipmi_ch = nothing 
    
    # Warmup code if we are doing a multi-sample timing run
    if nsamples > 1
        @info "Doing initial warm-up run"
        for event in events
            _ = inclusive_jets(jet_reconstruct(event, R = radius, p = p,
            strategy = strategy); ptmin = ptmin)
        end
    end
    
    # Threading?
    threads = Threads.nthreads()
    if (threads < 2) && ipmi 
        @error "Need at least 2 threads to perform parallel ipmi measurements of power"
        return 0.0
    end
    
    # Now setup timers and run the loop
    cummulative_time = 0.0
    cummulative_time2 = 0.0
    lowest_energy = 0.0
    ipmi_nrg = missing
    rapl_nrg = missing
    lowest_time = typemax(Float64)
    finaljets = Vector{Vector{PseudoJet}}(undef, threads)
    fj = Vector{Vector{FinalJet}}(undef, threads)
    
    if ipmi
        ipmi_ch = ipmi_chan(1)
	put!(ipmi_ch,0.0) 
    end
    
    for irun in 1:nsamples
        
        t_start = time_ns()
        #this was threaded in Graeme's version but can't be now because we want to do parallel ipmi measurement
        for event_counter ∈ 1:n_events * repeats
            event_idx = mod1(event_counter, n_events)
            my_t = Threads.threadid()
            inclusive_jets(jet_reconstruct(events[event_idx], R = radius, p = p,
            strategy = strategy), ptmin = ptmin)
        end

        t_stop = time_ns()

        dt_μs = convert(Float64, t_stop - t_start) * 1.e-3
        if nsamples > 1
            @info "$(irun)/$(nsamples) $(dt_μs)"
        end
        cummulative_time += dt_μs
        cummulative_time2 += dt_μs^2
        if dt_μs < lowest_time
            lowest_time = dt_μs
        end
    end

    if ipmi
	put!(ipmi_ch, 0.0)
	ipmi_nrg = take!(ipmi_ch)
        rapl_nrg = take!(ipmi_ch)
    end
    
    mean = cummulative_time / nsamples
    cummulative_time2 /= nsamples
    if nsamples > 1
        sigma = sqrt(nsamples / (nsamples - 1) * (cummulative_time2 - mean^2))
    else
        sigma = 0.0
    end
    mean /= n_events * repeats
    sigma /= n_events * repeats
    lowest_time /= n_events * repeats
    # Why also record the lowest time? 
    # 
    # The argument is that on a "busy" machine, the run time of an application is
    # always TrueRunTime+Overheads, where Overheads is a nuisance parameter that
    # adds jitter, depending on the other things the machine is doing. Therefore
    # the minimum value is (a) more stable and (b) reflects better the intrinsic
    # code performance.
    (lowest_time, impi_nrg/cummulative_time, rapl_nrg/cummulative_time)
end

function fastjet_jet_process_avg_time(input_file::AbstractString;
    ptmin::Float64 = 5.0,
    radius::Float64 = 0.4,
    p::Union{Real, Nothing} = nothing,
    algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
    strategy::RecoStrategy.Strategy,
    nsamples::Integer = 1,
    repeats::Integer = 1, #added to keep positional ordering same
    ipmi::Bool = false )

    
    # FastJet reader cannot handle gzipped files
    if endswith(input_file, ".gz")
        input_file = hepmc3gunzip(input_file)
    end
    
    # Set consistent algorithm and power
    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
    algorithm = algorithm)
    
    # @warn "FastJet timing not implemented yet"
    fj_bin = joinpath(@__DIR__, "..", "fastjet", "build", "fastjet-finder")
    fj_args = String[]
    push!(fj_args, "-p", string(p))
    push!(fj_args, "-s", string(strategy))
    push!(fj_args, "-R", string(radius))
    push!(fj_args, "--ptmin", string(ptmin))
    
    push!(fj_args, "-n", string(nsamples))
    @info "Fastjet command: $fj_bin $fj_args $input_file"
    
    t_start = UInt64(0)
    t_us = missing
    ipmi_nrg = 0.0
    ratp_nrg = 0.0

    ipmi_ch = nothing
    if ipmi
        t_start = time_ns()
        ipmi_ch = ipmi_chan(1) 
        put!(ipmi_ch, 0.0)
    end
 
    fj_output = read(`$fj_bin $fj_args $input_file`, String)

    if ipmi
	put!(ipmi_ch, 0.0) 
    	t_end = time_ns()
        t_us = (t_end - t_start) * 1e-3
    	ipmi_nrg = take!(ipmi_ch) 
    	ratp_nrg = take!(ipmi_ch)     
    end

    min = tryparse(Float64, match(r"Lowest time per event ([\d\.]+) us", fj_output)[1])
    if isnothing(min)
        @error "Failed to parse output from FastJet script"
        return (0.0,0.0,0.0)
    end
    (min, ipmi_nrg / t_us, ratp_nrg / t_us)
end

function python_jet_process_avg_time(backend::Backends.Code,
    input_file::AbstractString;
    ptmin::Float64 = 5.0,
    radius::Float64 = 0.4,
    p::Union{Real, Nothing} = nothing,
    algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
    strategy::RecoStrategy.Strategy,
    nsamples::Integer = 1,
    repeats::Integer = 1, #added to keep positional ordering same
    ipmi::Bool = false )

    
    # Python reader cannot handle gzipped files
    if endswith(input_file, ".gz")
        input_file = hepmc3gunzip(input_file)
    end
    
    # Set consistent algorithm and power
    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
    algorithm = algorithm)

    # There are some limitations in the Python code - only AntiKt is supported,
    # and the strategy has to be manually set
    if algorithm != JetAlgorithm.AntiKt
        @error "Only AntiKt is supported for Python backend"
        return 0.0
    end
    if strategy == RecoStrategy.Best
        @error "Best strategy is not supported for Python backend - use N2Plain or N2Tiled"
        return 0.0
    end

    py_args = String[]
    py_script = joinpath(@__DIR__, "..", "antikt-python", "src")
    if strategy == RecoStrategy.N2Plain
        py_script = joinpath(py_script, "antikt-basic.py")
    elseif strategy == RecoStrategy.N2Tiled 
        py_script = joinpath(py_script, "antikt-tiledN2cluster.py")
    else
        @error "Invalid strategy for Python: $strategy"
        return 0.0
    end

    # Accelerated or not?
    if backend == Backends.AkTNumPy
        py_args = ["--numba"]
    end
    
    # Ensure that all events are processed
    push!(py_args, "--maxevents", "100")

    push!(py_args, "--radius", string(radius))
    push!(py_args, "--ptmin", string(ptmin))
    
    push!(py_args, "--trials", string(nsamples))
    @info "Python command: $py_script $py_args $input_file"

    t_start = UInt64(0)
    t_us = missing
    ipmi_nrg = 0.0
    ratp_nrg = 0.0

    ipmi_ch = nothing
    if ipmi
        t_start = time_ns()
        ipmi_ch = ipmi_chan(1) 
        put!(ipmi_ch, 0.0)
    end
 
    py_output = read(`$py_script $py_args $input_file`, String)
   

    if ipmi
	put!(ipmi_ch, 0.0)
    	t_end = time_ns()
        t_us = (t_end - t_start) * 1e-3
    	ipmi_nrg = take!(ipmi_ch) 
    	ratp_nrg = take!(ipmi_ch)     
    end

 
    min = tryparse(Float64, match(r"Minimum time per event ([\d\.]+) us", py_output)[1])
    if isnothing(min)
        @error "Failed to parse output from Python script"
        return (0.0,0.0)
    end
    (min, ipmi_nrg/t_us, ratp_nrg/t_us)
end

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
    "--ptmin"
    help = "Minimum p_t for final jets (GeV)"
    arg_type = Float64
    default = 5.0
    
    "--radius", "-R"
    help = "Radius parameter for jet merging"
    arg_type = Float64
    default = 0.4
    
    "--algorithm", "-A"
    help = """Algorithm to use for jet reconstruction: $(join(JetReconstruction.AllJetRecoAlgorithms, ", "))"""
    arg_type = JetAlgorithm.Algorithm
    default = JetAlgorithm.AntiKt
    
    "--power", "-p"
    help = """Power value for jet reconstruction"""
    arg_type = Float64
    
    "--strategy", "-S"
    help = """Strategy for the algorithm, valid values: $(join(JetReconstruction.AllJetRecoStrategies, ", "))"""
    arg_type = RecoStrategy.Strategy
    default = RecoStrategy.Best
    
    "--nsamples", "-m"
    help = "Number of measurement points to acquire."
    arg_type = Int
    default = 16
    
    "--nsamples-override"
    help = "Override for sample number for specific event files"
    nargs = '+'
    
    "--repeats"
    help = "Run over whole event sample this number of times"
    arg_type = Int
    default = 1
    
    "--code"
    help = """Code backend to use for the jet reconstruction: $(join(AllCodeBackends, ", "))"""
    arg_type = Backends.Code
    default = Backends.JetReconstruction
    
    "--code-version"
    help = "Specific version string for the code backend used"
    arg_type = String
    default = "unknown"

    "--backend"
    help = """Backend used by the code - will be set automatically for Python and Julia, but Fastjet may benefit from being set manually (e.g., gcc or clang)"""
    arg_type = String

    "--backend-version"
    help = """Specific version string for the backend used - will be set automatically for Python and Julia, but Fastjet may benefit from being set manually"""
    arg_type = String

    "--info"
    help = "Print info level log messages"
    action = :store_true
    
    "--debug"
    help = "Print debug level log messages"
    action = :store_true
    
    "--results"
    help = """Write results in CSV format to this directory/file. If a directory is given, a file named 'BACKEND-ALGORITHM-STRATEGY-RADIUS.csv' will be created."""
    
    "--ipmi"
    help = "Attempt to measure power from IPMI (system) and RATP (CPU)"
    action = :store_true

    "files"
    help = "HepMC3 event files in to process or CSV file listing event files"
    required = true
    nargs = '+'
end
return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    if args[:debug]
        logger = ConsoleLogger(stdout, Logging.Debug)
    elseif args[:info]
        logger = ConsoleLogger(stdout, Logging.Info)
    else
        logger = ConsoleLogger(stdout, Logging.Warn)
    end
    global_logger(logger)
    
    # If we have a CSV file, open it and read
    if endswith(args[:files][1], ".csv")
        if length(args[:files]) > 1
            println("When CSV input is given, only one file can be used")
            exit(1)
        end
        hepmc3_files_df = CSV.read(args[:files][1], DataFrame)
        input_file_path = dirname(args[:files][1])
        input_files = hepmc3_files_df[:, :File]
        input_files_full = map(fname -> joinpath(input_file_path, fname), input_files)
        hepmc3_files_df[:, :File_path] = input_files_full
    else
        # It's just a plain list of files
        input_files_full = args[:files]
        hepmc3_files_df = DataFrame("File_path" => input_files_full)
        hepmc3_files_df[:, :File] = map(basename, hepmc3_files_df[:, :File_path])
        hepmc3_files_df[:, :mean_particles] .= -1
    end

    # Get consistent algorithm and power here, so that missing values are filled
    (power, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = args[:power], algorithm = args[:algorithm])
    
    event_timing = Float64[]
    ipmi_powering = Union{Missing,Float64}[]
    ratp_powering = Union{Missing,Float64}[]

    n_samples = Int[]
    for event_file in hepmc3_files_df[:, :File_path]
        if event_file in args[:nsamples_override]
            samples = tryparse(Int, args[:nsamples_override][1])
            @info "Overriding samples for $(event_file) to $samples"
        else
            samples = args[:nsamples]
        end
        push!(n_samples, samples)
        
        if args[:code] == Backends.JetReconstruction
            # Try to read events into the correct type!
            if JetReconstruction.is_ee(args[:algorithm])
                JetType = EEJet
            else
                JetType = PseudoJet
            end
            events::Vector{Vector{JetType}} = read_final_state_particles(event_file; T = JetType)
            time_per_event,ipmi_power,ratp_power = julia_jet_process_avg_time(events; ptmin = args[:ptmin],
            radius = args[:radius],
            algorithm = args[:algorithm],
            p = args[:power],
            strategy = args[:strategy],
            nsamples = samples, repeats = args[:repeats], ipmi=args[:ipmi])
        elseif args[:code] == Backends.Fastjet
            time_per_event,ipmi_power,ratp_power = fastjet_jet_process_avg_time(event_file; ptmin = args[:ptmin],
            radius = args[:radius],
            algorithm = args[:algorithm],
            p = args[:power],
            strategy = args[:strategy],
            nsamples = samples, ipmi=args[:ipmi])
        elseif args[:code] in (Backends.AkTPython, Backends.AkTNumPy)
            time_per_event,ipmi_power,ratp_power = python_jet_process_avg_time(args[:code], event_file; ptmin = args[:ptmin],
            radius = args[:radius],
            algorithm = args[:algorithm],
            p = args[:power],
            strategy = args[:strategy],
            nsamples = samples, ipmi=args[:ipmi])
        end
        
        push!(event_timing, time_per_event)
        push!(ipmi_powering, ipmi_power)
        push!(ratp_powering, ratp_power)
    end
    # Add results to the DataFrame
    hepmc3_files_df[:, :n_samples] = n_samples
    hepmc3_files_df[:, :time_per_event] = event_timing
    hepmc3_files_df[:, :ipmi_power] = ipmi_powering
    hepmcs_files_df[:, :ratp_power] = ratp_powering

    # Decorate the DataFrame with the metadata of the runs
    hepmc3_files_df[:, :code] .= args[:code]
    hepmc3_files_df[:, :code_version] .= args[:code_version]
    hepmc3_files_df[:, :algorithm] .= algorithm
    hepmc3_files_df[:, :strategy] .= args[:strategy]
    hepmc3_files_df[:, :radius] .= args[:radius]
    hepmc3_files_df[:, :power] .= power

    backend = isnothing(args[:backend]) ? determine_backend(args[:code]) : args[:backend]
    backend_version = isnothing(args[:backend_version]) ? determine_backend_version(args[:code]) : args[:backend_version]
    hepmc3_files_df[:, :backend] .= backend
    hepmc3_files_df[:, :backend_version] .= backend_version

    println(hepmc3_files_df)
    
    # Write out the results
    if !isnothing(args[:results])
        if isdir(args[:results])
            results_file = joinpath(args[:results],
            "$(args[:code])_$(args[:code_version])_$(args[:algorithm])_" *
            "$(args[:strategy])_R$(args[:radius])_P$(power)_$(backend)_$(backend_version).csv")
        else
            results_file = args[:results]
        end
        @info "Writing results to $(results_file)"
        CSV.write(results_file, hepmc3_files_df)
    end
    
    nothing
end

main()
