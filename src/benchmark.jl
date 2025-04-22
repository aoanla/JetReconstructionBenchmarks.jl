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
    output = read(`ipmitool dcmi power reading`, String)
    t = time_ns()
    s = tryparse(Float64, match(r"Instantaneous power reading:\s+([\d\.]+) Watts", output)[1])
    (s, t)
end

function ipmi_chan(nsamples:Integer)
    Channel{Float64}(spawn=true) do ch 
        sample_time = [(0.0,0.0), (0.0,0.0)] 
        idx = 1
        energy_nj = 0.0

        rdy() = (isready(ch))
    
        for irun in 1:nsamples
            take!(ch) #block until we tell it to start
            sample_time[idx] = ipmi_command()
            idx = 3 - idx #toggle 1 and 2
            while !rdy()
                    timedwait(rdy, 2) #sample every 2 seconds or until we're told to stop
                    sample_time[idx] = ipmi_command()
                    idx = 3 - idx
                    #trapezium rule
                    energy_nj += (sample_time[1][1]+sample_time[2][1])*abs(sample_time[2][1]-sample_time[2][2])/2.0
            end
            take!(ch) #take the value passed to us to stop our counter 
            put!(energy_nj * 1.e-9 ) #Joules
            energy_nj = 0.0
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
    nrg = 0.0
    lowest_time = typemax(Float64)
    finaljets = Vector{Vector{PseudoJet}}(undef, threads)
    fj = Vector{Vector{FinalJet}}(undef, threads)
    
    if ipmi
        ipmi_ch = ipmi_chan(nsamples) 
    end
    
    for irun in 1:nsamples
        
        if ipmi    
            put!(ipmi_ch,0.0)
        end
        
        t_start = time_ns()
        #this was threaded in Graeme's version but can't be now because we want to do parallel ipmi measurement
        for event_counter ∈ 1:n_events * repeats
            event_idx = mod1(event_counter, n_events)
            my_t = Threads.threadid()
            inclusive_jets(jet_reconstruct(events[event_idx], R = radius, p = p,
            strategy = strategy), ptmin = ptmin)
        end

        if ipmi
            put!(ipmi_ch,0.0)
        end
        
        t_stop = time_ns()

        if ipmi
            nrg = take!(ipmi_ch)
        end
        
        dt_μs = convert(Float64, t_stop - t_start) * 1.e-3
        if nsamples > 1
            @info "$(irun)/$(nsamples) $(dt_μs)"
        end
        cummulative_time += dt_μs
        cummulative_time2 += dt_μs^2
        if dt_μs < lowest_time
            lowest_time = dt_μs
            if ipmi
                lowest_energy = nrg
            end
        end
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
    (lowest_time, lowest_energy)
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

    if ipmi && (nsamples > 1)
        @error "Cannot perform ipmi power sampling for multiple samples for fastjet"
        return 0.0
    end
    
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
    
    ipmi_ch = nothing
    if ipmi
        ipmi_ch = ipmi_chan(nsamples) 
        put!(ipmi_ch, 0.0)
    end 
    fj_output = read(`$fj_bin $fj_args $input_file`, String)
    put!(ipmi_ch, 0.0)
    energy = take!(ipmi_ch)
    
    min = tryparse(Float64, match(r"Lowest time per event ([\d\.]+) us", fj_output)[1])
    if isnothing(min)
        @error "Failed to parse output from FastJet script"
        return (0.0,0.0)
    end
    (min, energy)
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

    if ipmi && (nsamples > 1)
        @error "Cannot perform ipmi power sampling for multiple samples for python"
        return 0.0
    end
    
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

    ipmi_ch = nothing
    if ipmi
        ipmi_ch = ipmi_chan(nsamples) 
        put!(ipmi_ch, 0.0)
    end 

    py_output = read(`$py_script $py_args $input_file`, String)
    put!(ipmi_ch, 0.0)
    energy = take!(ipmi_ch)
    
    min = tryparse(Float64, match(r"Minimum time per event ([\d\.]+) us", py_output)[1])
    if isnothing(min)
        @error "Failed to parse output from Python script"
        return (0.0,0.0)
    end
    (min, energy)
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
            time_per_event = julia_jet_process_avg_time(events; ptmin = args[:ptmin],
            radius = args[:radius],
            algorithm = args[:algorithm],
            p = args[:power],
            strategy = args[:strategy],
            nsamples = samples, repeats = args[:repeats])
        elseif args[:code] == Backends.Fastjet
            time_per_event = fastjet_jet_process_avg_time(event_file; ptmin = args[:ptmin],
            radius = args[:radius],
            algorithm = args[:algorithm],
            p = args[:power],
            strategy = args[:strategy],
            nsamples = samples)
        elseif args[:code] in (Backends.AkTPython, Backends.AkTNumPy)
            time_per_event = python_jet_process_avg_time(args[:code], event_file; ptmin = args[:ptmin],
            radius = args[:radius],
            algorithm = args[:algorithm],
            p = args[:power],
            strategy = args[:strategy],
            nsamples = samples)
        end
        
        push!(event_timing, time_per_event)
    end
    # Add results to the DataFrame
    hepmc3_files_df[:, :n_samples] = n_samples
    hepmc3_files_df[:, :time_per_event] = event_timing

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
