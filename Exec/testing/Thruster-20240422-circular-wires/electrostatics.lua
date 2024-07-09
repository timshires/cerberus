verbosity = 2
linear_solver_verbosity = 1

cfl = 0.5

-- === DEFINE PROBLEM ===

lightspeed = 1000.0
Larmor = 1.0
Debye = 1.0

function ion_number_density(dat)
    x = dat['x'] - 4.0
    y = dat['y']
    
    r = math.sqrt(x*x + y*y)
    
    if (r < 1.0) then
        return 1.0
    else
        return 1.0
    end
end

-- === DEFINE STATES ===

states = {

    field = {
        type='field',
        reconstruction='O6',
        static=1,
    },
    
    ion = {
        type='hydro',
        gas={
          type='thermally_perfect',
          mass=1.0,
          charge=1.0,
          gamma=5/3,
        },
        reconstruction='minmod', 
        flux='HLLC',
        refinement={name='hydro_gradient', rho=0.1},
        value = {
            nd = ion_number_density,
            p =   1.0,
        },
    },

    electron = {
        type='hydro',
        gas={
          type='thermally_perfect',
          mass=1.0,
          charge=-1.0,
          gamma=5/3,
        },
        reconstruction='minmod',  
        flux='HLLC',
        refinement={name='hydro_gradient', rho=0.1},
        value = {
            nd = 1.0,
            p =  1.0,
        },
    },
}

-- === SOURCE TERMS ===

actions = {
    hydro_fluxes = {
        type = 'CTU',
        corner_transport=true,
        states = {'ion', 'electron'},
    },

    plasma={
        type='Lorentz',
        states = {'ion', 'electron', 'field'},
     },

     fields = {
        type='elliptic',
        projection=0,
        state = 'field',
     },
}

-- === GEOMETRY ===

refine_cutcells = true


function make_circles(x,y,collection)

    local d, dd

    for i, v in ipairs(collection) do
      dd = v[1]^2 - ((x-v[2])^2 + (y-v[3])^2)
      if (i == 1) then
          d = dd
      end
      d = math.max(d, dd)
    end

    return -d

end

function cylinder(x,y)

    local circles = {
      {7.5, 0.0, 0.0, 0.0},
    }

    return make_circles(x,y,circles)
end

function anode(x,y)

    local circles = {
      {1.0, 0.0, 4.0, 0.0},
    }

    return -make_circles(x,y,circles)
end

function cathode(x,y)

    local circles = {
      {1.0, 0.0, -4.0, 0.0},
    }

    return -make_circles(x,y,circles)
end


embedded_boundaries = {

    shell = {
      geom=cylinder,
      bcs={
        field={types={'scalar_potential', 'vector_potential'}, phi=0.0, A1=0.0, A2=-1.0},
        ion={type='slip_wall'},
        electron={type='slip_wall'},
      },
      boolean_operation='and',
    },
    
    anode = {
      geom=anode,
      bcs={
        field={types={'scalar_potential', 'vector_potential'}, phi=-0.01, A1=0.0, A2=1.0},
        ion={type='slip_wall'},
        electron={type='slip_wall'},
      },
      boolean_operation='and',
    },
    
    cathode = {
      geom=cathode,
      bcs={
        field={type='scalar_potential', phi=-0.0},
        ion={type='slip_wall'},
        electron={type='slip_wall'},
        --electron={type='dirichlet',rho=1.0, x_vel=0.0},
      },
      boolean_operation='and',
    },

}

