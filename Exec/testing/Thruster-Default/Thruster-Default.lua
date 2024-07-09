-- ======== PROBLEM ==========
--# Assorted constants for calculations - dont change #--
-- mu_0_dim = 1.25663706e-6; 
-- ep_0_dim = 8.85418782e-12;
-- kb = 1.38064852e-23;
-- ref_lightspeed = 299792458.0; --dimensional speed of light 

--# Reference Values 
-- q_ref = 1.60217662e-19 -- Coulombs

-- m_p = 1.6726219000e-27 -- Proton mass
-- mfac = 1 --  1e17  --TODO problem specific for ionisation fraction or for stability? 
-- m_ref = mfac * m_p

-- n_ref = 1.0e23

-- REFERENCE VALUES 
-- ref_length = 0.1 -- m 
-- ref_mass = m_ref
-- ref_density = n_ref * ref_mass -- use the number density multipled by the ref mass 


--ref_temp = 273.0
-- ref_temp = 1.2097998979e+06

--Larmor = 0.1
--Debye = 0.1
-- Larmor = 1.0e-2 --TODO google what these are/represent --TODO HW
-- Debye  = 2.4e-6


-- MLMG solver friendly numbers --- use as a starting point and worry about 
--reference parameters once things are running ...

lightspeed = 1000.0
Larmor = 1.0
Debye = 1.0

-- phi_mag = 1 -- target D_x nondim
-- A_mag   = 0 -- target B_z non-dim

-- === SETTINGS ===
verbosity = 2
linear_solver_verbosity = 1
cfl             = 0.5 -- 0.25
-- do_face_sources = 0
-- do_CTU          = 1

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

-- === functions ===
--  Fluids initial condition functions (based on spatial distribution of cells)
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

-- function tracer(dat)
--   x = dat['x']
--   y = dat['y']

--   if (x < 2.) then
--     t = 0.;
--   else
--     t = 1.0;
--   end

--   return t
-- end 

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

-- === Define the states ===

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
    -- value = {
    --   x_D = 0.0,
    --   y_D = 0.0,
    --   z_D = 0.0,
    --   x_B = 0.0,
    --   y_B = 0.0,
    --   z_B = 0.0},
    -- bc={x={lo={fill_D_bc='symmetry',
    --            fill_B_bc='symmetry',
    --            fill_psi_bc='symmetry',
    --            fill_phi_bc='symmetry',
    --            fill_ep_bc='symmetry',
    --            fill_mu_bc='symmetry',
    --           },
    --        hi={fill_D_bc='symmetry',
    --            fill_B_bc='symmetry',
    --            fill_psi_bc='symmetry',
    --            fill_phi_bc='symmetry',
    --            fill_ep_bc='symmetry',
    --            fill_mu_bc='symmetry',
    --           },
    --        },
    --     y={lo={fill_D_bc='symmetry',
    --            fill_B_bc='symmetry',
    --            fill_psi_bc='symmetry',
    --            fill_phi_bc='symmetry',
    --            fill_ep_bc='symmetry',
    --            fill_mu_bc='symmetry',
    --           },
    --        hi={fill_D_bc='symmetry',
    --            fill_B_bc='symmetry',
    --            fill_psi_bc='symmetry',
    --            fill_phi_bc='symmetry',
    --            fill_ep_bc='symmetry',
    --            fill_mu_bc='symmetry',
    --           },
    --        },
    --   },
  },
}

-- === SOURCE TERMS ===
--TODO put back in first 
actions = {
  plasma = {
    type='Lorentz',
    -- solver   = 'implicit',
    -- options  = {order=3},
    states = {'ion', 'electron','field'}
  },

  hydro_fluxes = {
    type = 'CTU',
    corner_transport=true,
    states = {'ion', 'electron'},
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

--TODO put back in slowly once running 
-- collection = {{{x_lo, x_hi},{y_lo,y_hi}}, ... }
-- function make_rectangles(x,y,collection)

--   local d, d1, dx, dy, dx_, dy_

--   local coords = {x,y}

--   for i, v in ipairs(collection) do
--     dx = math.max(coords[1] - v[1][2], v[1][1] - coords[1])
--     dy = math.max(coords[2] - v[2][2], v[2][1] - coords[2])

--     dx_ = math.max(dx, 0.0)
--     dy_ = math.max(dy, 0.0)

--     d1 = math.sqrt(dx_*dx_ + dy_*dy_) + math.min(0.0, math.max(dx, dy))

--     if (i == 1) then
--       d = d1
--     else 
--       d = math.min(d, d1)
--     end
--   end

--   return d    
-- end

-- -- {xlo, xhi}, {ylo, yhi}
-- x_grid = {2,2.1}

-- --  # establish end grid for ion engine, formed basically using block 1 and 2.
-- -- function block1(x,y)
-- --   local rect = {
-- --     {x_grid,{0.4, 0.8}},
-- --   }
-- --   return make_rectangles(x,y,rect)
-- -- end
-- -- function block2(x,y)
-- --   local rect = {
-- --     {x_grid,{1.2, 1.6}},
-- --   }
-- --   return make_rectangles(x,y,rect)
-- -- end

-- function block3(x,y)
--   local rect = {
--     {{3,3.1},{0, 1}},
--   }
--   return make_rectangles(x,y,rect)
-- end
-- function block4(x,y)
--   local rect = {
--     {{1,1.1},{0, 1}},
--   }
--   return make_rectangles(x,y,rect)
-- end

function phi_f(dat)
  local phival
  x = dat['x']
  y = dat['y']

  phival=-x*phi_mag
  return phival
end

embedded_boundaries = {
  -- -- block1 = {
  -- --   geom = block1,
  -- --   bcs={
  -- --     field={
  -- --       types={'scalar_potential'}, 
  -- --       phi=phi_mag
  -- --     },
  -- --     -- ion={
  -- --     --   type='slip_wall' -- try omitting to see for non-interaction
  -- --     -- },
  -- --     -- electron={
  -- --     --   type='slip_wall'
  -- --     -- },
  -- --   },
  -- --   boolean_operation='and',
  -- --   inside=0,
  -- -- },
  -- -- block2 = {
  -- --   geom = block2,
  -- --   bcs={
  -- --     field={
  -- --       types={'scalar_potential'}, 
  -- --       phi=phi_mag
  -- --     },
  -- --     ion={
  -- --       type='slip_wall'
  -- --     },
  -- --     electron={
  -- --       type='slip_wall'
  -- --     },
  -- --   },
  -- --   boolean_operation='and',
  -- --   inside=0,
  -- -- },
  -- block3 = {
  --   geom = block3,
  --   bcs={
  --     field={
  --       types={'scalar_potential'}, 
  --       phi=-0.01
  --     },
  --     -- ion={
  --     --   type='slip_wall'
  --     -- },
  --     -- electron={
  --     --   type='slip_wall'
  --     -- },
  --   },
  --   boolean_operation='and',
  --   -- inside=0, -- set to one to matchdefaut used in EM-electrostatic??
  -- },
  -- block4 = {
  --   geom = block4,
  --   bcs={
  --     field={
  --       types={'scalar_potential'}, 
  --       phi=0.01
  --     },
  --     -- ion={
  --     --   type='slip_wall'
  --     -- },
  --     -- electron={
  --     --   type='slip_wall'
  --     -- },
  --   },
  --   boolean_operation='and',
  --   -- inside=0,
  -- },
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
--TODO === EB put back in when deubg done ===

-- === PLOTTING ===

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

