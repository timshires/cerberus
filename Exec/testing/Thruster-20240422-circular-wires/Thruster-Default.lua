-- ======== PROBLEM ==========
-- MLMG solver friendly numbers --- use as a starting point and worry about 
--reference parameters once things are running ...
lightspeed = 1000.0
Larmor = 1.0
Debye = 1.0

phi_mag = 1 -- target D_x nondim
A_mag   = 0 -- target B_z non-dim

-- === SETTINGS ===
verbosity       = 2
linear_solver_verbosity = 1
cfl             = 0.5 -- 0.25
-- do_face_sources = 0
-- do_CTU          = 1

-- === DEFINE STATES ===
density       = 1.0
mass_ion      = 1
mass_electron = 1
gam           = 5./3.
pressure      = 1

-- computations
-- p0   = pressure
-- rho0 = density
-- a0   = math.sqrt(gam*p0/rho0)
-- u0   = 0.0

-- override:
-- p1   = 1. --4.13291e-06
-- rho1 = 1. --0.014
-- u1   = 0.01

-- === functions ===
--  Fluids initial condition functions (based on spatial distribution of cells)
function density(dat)

    x = dat['x']
    y = dat['y']

    rho = rho1

    return rho
end


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


-- function ion_density(dat)
--     return density(dat)
-- end

-- function electron_density(dat)
--     return mass_electron*density(dat)/mass_ion -- we are balancing number density - think about this and explain why HW - (relative density to ion?)
-- end


-- function pressure(dat)
--   x = dat['x']
--   y = dat['y']
--   return p1
-- end

-- function velocity_x(dat)
--   x = dat['x']
--   y = dat['y']
--   return u1
-- end

-- === Define the states ===

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


--embedded boundary definition 

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

-- === PLOTTING ===

--[[
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
--]]
