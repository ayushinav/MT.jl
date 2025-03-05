params_HZK2011 = (
    mechs = (
        diff=(
            A=10.0f0^76.0f-1,  # Preexponential for coble diffusion creep
            Q=375.0f3,   # Activation energy for coble diffusion creep
            V=10.0f-6,   # Activation volume for coble diffusion creep
            p=3,       # Grain size exponent
            alf=25,    # Melt factor
            r=0,       # Water fugacity exponent
            n=1,       # Stress exponent
            ϕ_c=1.0f-5,
            x_ϕ_c=5
        ), 
        
        disl=(
            A=1.1f5,   # Preexponential
            Q=530.0f3,   # Activation energy
            V=15.0f-6,   # Activation volume
            n=35.0f-1,     # Stress exponent
            p=0,       # Grain size exponent
            alf=30,    # Melt factor
            r=0,       # Water fugacity exponent
            ϕ_c=1.0f-5,
            x_ϕ_c=1
        ), 
        
        gbs=(
            A=10.0f0^4.8f0,  # Preexponential for GBS-disl creep
            Q=445.0f3,   # Activation energy for GBS-disl creep
            V=15.0f-6,   # Activation volume
            p=73.0f-2,    # Grain size exponent
            n=29.0f-1,     # Stress exponent
            alf=35,    # Melt factor
            r=0,       # Water fugacity exponent
            ϕ_c=1.0f-5,
            x_ϕ_c=25.0f-1
        ),
        ),
    p_dep_calc = true,
    melt_enhancement = false
    
    )