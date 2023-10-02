module QECC_Tests

export plot_code_performance, evaluate_degen_code_decoder, evaluate_code_decoder, evaluate_code_decoder_noisy_circuit

 using QuantumClifford, CairoMakie
 using QuantumClifford.ECC: naive_syndrome_circuit

 function plot_code_performance(error_rates, post_ec_error_rates; title="")
    f = Figure(resolution=(500,300))
    ax = f[1,1] = Axis(f, xlabel="single (qu)bit error rate",title=title)
    ax.aspect = DataAspect()
    lim = max(error_rates[end],post_ec_error_rates[end])
    lines!([0,lim], [0,lim], label="single bit", color=:black)
    plot!(error_rates, post_ec_error_rates, label="after decoding", color=:black)
    xlims!(0,lim)
    ylims!(0,lim)
    f[1,2] = Legend(f, ax, "Error Rates")
    f
 end;

 function evaluate_code_decoder(H,lookup_table,p; samples=10_000)
    constraints, bits = size(H)
    decoded = 0 # Counts correct decodings
    for sample in 1:samples
        error = rand(bits) .< p                    # Generate random error
        syndrome = (H*error) .% 2                  # Apply that error to your physical system and get syndrome
        guess = get(lookup_table,syndrome,nothing) # Decode the syndrome
        if guess==error                            # Check if you were right
            decoded += 1
        end
    end
    1 - decoded / samples
 end;

 function evaluate_code_decoder(code::Stabilizer,lookup_table,p; samples=10_000)
    constraints, qubits = size(code)
    decoded = 0 # Counts correct decodings
    for sample in 1:samples
        # Generate random error
        error = random_pauli(qubits,p/3,nophase=true)
        # Apply that error to your physical system
        # and get syndrome
        syndrome = comm(error, code)
        # Decode the syndrome
        guess = get(lookup_table,syndrome,nothing)
        # check if you were right
        if guess==error
            decoded += 1
        end
    end
    1 - decoded / samples
 end;

 function evaluate_degen_code_decoder(code::Stabilizer,lookup_table,p; samples=10_000)
    constraints, qubits = size(code)
    full_tableau = MixedDestabilizer(code)
    logicals = vcat(logicalxview(full_tableau),logicalzview(full_tableau))
    decoded = 0 # Counts correct decodings
    for sample in 1:samples
        # Generate random error
        error = random_pauli(qubits,p/3,nophase=true)
        # Apply that error to your physical system
        # and get syndrome
        syndrome = comm(error, code)
        # Decode the syndrome
        guess = get(lookup_table,syndrome,nothing)
        # Check if the suggested error correction
        # corrects the error or if it is equivalent
        # to a logical operation
        if !isnothing(guess) && all(==(0x0), comm(guess*error, code)) && all(==(0x0), comm(guess*error, logicals))
            decoded += 1
        end
    end
    1 - decoded / samples
 end; 

 function evaluate_code_decoder_noisy_circuit(code::Stabilizer,lookup_table,p,q; samples=10_000)
    s, n = size(code)
    constraints, qubits = size(code)
    initial_state = state = Register(MixedDestabilizer(code ⊗ one(Stabilizer, s)),zeros(Bool,6))
    initial_state.stab.rank = qubits+constraints # TODO hackish and ugly, needs fixing
    no_error_state = canonicalize!(stabilizerview(traceout!(copy(initial_state),n+1:n+constraints))) # TODO rather ugly
    # prepare measurement circuit
    syndrome_circuit = naive_syndrome_circuit(code);
    noisy_circuit = make_noisy(syndrome_circuit, UnbiasedUncorrelatedNoise(q/3))
    mem_noise_op = NoiseOp(UnbiasedUncorrelatedNoise(p/3), 1:qubits)
    pad_ancilla_qubits = zero(PauliOperator, constraints) # TODO rather ugly
    decoded = 0 # Counts correct decodings
    for sample in 1:samples
        state = copy(initial_state)
        # Apply the initial memory error to your physical system
        apply!(state, mem_noise_op)
        # Run the syndrome measurement circuit
        mctrajectory!(state, noisy_circuit)
        syndrome = UInt8.(bitview(state))
        # Decode the syndrome
        guess = get(lookup_table,syndrome,nothing)
        if isnothing(guess)
            continue
        end
        # Apply the suggested correction
        apply!(state, guess ⊗ pad_ancilla_qubits)
        # Check for errors
        if no_error_state == canonicalize!(stabilizerview(traceout!(state,n+1:n+constraints))) # TODO rather ugly
            decoded += 1
        end
    end
    1 - decoded / samples
 end;

end