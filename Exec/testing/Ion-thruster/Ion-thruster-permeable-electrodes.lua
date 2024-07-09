verbosity = 2
linear_solver_verbosity = 1

cfl = 0.5

-- === DEFINE PROBLEM ===

lightspeed = 1000.0
Larmor = 1.0
Debye = 1.0


-- === DEFINE STATES ===
density       = 1.0
mass_ion      = 1 -- in electrostatic it is set to 1
mass_electron = 1 -- in electrostatic it is set to 1, same as electron
gam           = 5./3.
pressure      = 0.5

-- computations
p0   = pressure
rho0 = density
a0   = math.sqrt(gam*p0/rho0)
u0   = 0.0

-- override:
p1   = 1. --4.13291e-06
rho1 = 1. --0.014
u1   = 0.01

function density(dat)

  x = dat['x']
  y = dat['y']

  rho = rho1

  return rho
end

function ion_density(dat)
  return density(dat)
end

function electron_density(dat)
  return mass_electron*density(dat)/mass_ion -- we are balancing number density - think about this and explain why HW - (relative density to ion?)
end

function pressure(dat)
  x = dat['x']
  y = dat['y']
  return p1
end

function velocity_x(dat)
  x = dat['x']
  y = dat['y']
  return u1
end

-- === DEFINE STATES ===

states = {
  electron = {
    type = 'hydro',
    gas = {
      type   = 'thermally_perfect', 
      mass   = mass_electron,
      charge = -1.0,
      gamma  = 5./3.,
    },
    reconstruction        = 'minmod',
    flux                  = 'HLLC',
    refinement={name='hydro_gradient', rho = 0.1},
    value = {
      rho   = electron_density,
      x_vel = velocity_x,
      p     = pressure,
      alpha = tracer,
    },
    bc = {
      x={
        lo={
            gamma  = 5./3.,
            cp = 100,
            fill_hydro_bc = 'inflow',
            rho=electron_density, x_vel=velocity_x, y_vel=0.0, z_vel=0.0, T=0.0041329148800799, p=pressure, alpha=tracer,
        },
      },
    },
  },
  ion = {
    type ='hydro',
    gas  = {
      type   = 'thermally_perfect',
      mass   = mass_ion,
      charge = 1.0,
      gamma  = 5/3,
    },
    reconstruction = 'minmod', 
    flux           = 'HLLC',
    refinement     = {name='hydro_gradient', rho = 0.1},
    value = {
      rho   = ion_density,
      x_vel = velocity_x,
      p     = pressure,
      alpha = tracer,
    },
    bc = {
      x={
        lo={
          gamma  = 5./3.,
          cp = 100,
            fill_hydro_bc = 'inflow',
            rho=ion_density, x_vel=velocity_x, y_vel=0.0, z_vel=0.0, T=0.0041329148800799, p=pressure, alpha=tracer,
        },
      },
    },
  },
  field = {
    type           = 'field',
    static         = 1,
    reconstruction = 'O6',
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


function make_rectangles(x,y,collection)

  local d, d1, dx, dy, dx_, dy_

  local coords = {x,y}

  for i, v in ipairs(collection) do
    dx = math.max(coords[1] - v[1][2], v[1][1] - coords[1])
    dy = math.max(coords[2] - v[2][2], v[2][1] - coords[2])

    dx_ = math.max(dx, 0.0)
    dy_ = math.max(dy, 0.0)

    d1 = math.sqrt(dx_*dx_ + dy_*dy_) + math.min(0.0, math.max(dx, dy))

    if (i == 1) then
      d = d1
    else 
      d = math.min(d, d1)
    end
  end

  return d    
end

-- {xlo, xhi}, {ylo, yhi}
x_grid = {2,2.1}

--  establish end grid for ion engine, formed basically using block 1 and 2.
function block1(x,y)
  local rect = {
    {x_grid,{0.4, 0.8}},
  }
  return make_rectangles(x,y,rect)
end
function block2(x,y)
  local rect = {
    {x_grid,{1.2, 1.6}},
  }
  return make_rectangles(x,y,rect)
end

function block3(x,y)
  local rect = {
    {{3,3.1},{0.5, 1}},
  }
  return -make_rectangles(x,y,rect) -- adding minus sign made it run...?
end
function block4(x,y)
  local rect = {
    {{1,1.1},{0.5, 1}},
  }
  return -make_rectangles(x,y,rect) -- adding minus sign made it run...?
end


embedded_boundaries = {

    anode = {
      geom=block3,
      bcs={
        field={types={'scalar_potential', 'vector_potential'}, phi=-0.01, A1=0.0, A2=1.0},
        -- ion={type='slip_wall'},
        -- electron={type='slip_wall'},
      },
      boolean_operation='and',
    },
    
    cathode = {
      geom=block4,
      bcs={
        field={type='scalar_potential', phi=-0.0},
        -- ion={type='slip_wall'},
        -- electron={type='slip_wall'},
        --electron={type='dirichlet',rho=1.0, x_vel=0.0},
      },
      boolean_operation='and',
    },

}


function outside(dat)
  if ((dat['vfrac-electron'] == 0.0) or (dat['vfrac-ion'] == 0.0)) then
      return true
  else 
      return false
  end
end

function charge_density(dat)

  if (outside(dat)) then
      return 0.0
  end

  local cd_e = dat['charge-electron']*dat['rho-electron']/dat['mass-electron']
  local cd_i = dat['charge-ion']*dat['rho-ion']/dat['mass-ion']
  
  return (cd_e + cd_i)*(dat['Larmor']/dat['Debye']^2)
end

function div_D(dat)

  if (outside(dat)) then
      return 0.0
  end
  
  return dat['x_D-field-dx'] + dat['y_D-field-dy']
end

function err(dat)
  local cd = charge_density(dat)
  local div = div_D(dat)
  return cd - div
end

plot = {
  variables = {
      'all',
      'vfrac-electron',
      'vfrac-ion',
      'charge-electron',
      'charge-ion',
      'rho-electron',
      'rho-ion',
      'mass-electron',
      'mass-ion',
      'x_D-field-dx',
      'y_D-field-dy',
      'alpha_0-electron',
  },
  functions = {
      charge_density=charge_density,
      div_D=div_D,
      err_div_D=err,
  },
}

