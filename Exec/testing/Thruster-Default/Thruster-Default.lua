-- ======== PROBLEM ==========
--# Assorted constants for calculations - dont change #--
mu_0_dim = 1.25663706e-6; 
ep_0_dim = 8.85418782e-12;
kb = 1.38064852e-23;
ref_lightspeed = 299792458.0; --dimensional speed of light 

--# Reference Values 
q_ref = 1.60217662e-19 -- Coulombs

m_p = 1.6726219000e-27 -- Proton mass
mfac = 1. --  1e17  --TODO problem specific for ionisation fraction or for stability? 
m_ref = mfac * m_p

n_ref = 1.0e23

-- REFERENCE VALUES 
ref_length = 0.1 -- m 
ref_mass = m_ref
ref_density = n_ref * ref_mass -- use the number density multipled by the ref mass 


--ref_temp = 273.0
ref_temp = 1.2097998979e+06

--Larmor = 0.1
--Debye = 0.1
Larmor = 1.0e-2 --TODO google what these are/represent --TODO HW
Debye  = 2.4e-6

phi_mag = -0.001 -- target D_x nondim
A_mag   = 0 -- target B_z non-dim

-- === SETTINGS ===
verbosity       = 1
cfl             = 0.25
do_face_sources = 0
do_CTU          = 1

-- === DEFINE STATES ===
density       = 1.0
mass_ion      = 1
mass_electron = 0.01
gam           = 5./3.
pressure      = 0.5

-- computations
p0   = pressure
rho0 = density
a0   = math.sqrt(gam*p0/rho0)
u0   = 0.0

-- override:
p1   = 4.13291e-06
rho1 = 0.014
u1   = 1

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

function tracer(dat)
  x = dat['x']
  y = dat['y']

  if (x < 2.) then
    t = 0.;
  else
    t = 1.0;
  end

  return t
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
    value = {
      x_D = phi_mag,
      y_D = 0.0,
      z_D = 0.0,
      x_B = 0.0,
      y_B = 0.0,
      z_B = 0.0},
    bc={x={lo={fill_D_bc='symmetry',
               fill_B_bc='symmetry',
               fill_psi_bc='symmetry',
               fill_phi_bc='symmetry',
               fill_ep_bc='symmetry',
               fill_mu_bc='symmetry',
              },
           hi={fill_D_bc='symmetry',
               fill_B_bc='symmetry',
               fill_psi_bc='symmetry',
               fill_phi_bc='symmetry',
               fill_ep_bc='symmetry',
               fill_mu_bc='symmetry',
              },
           },
        y={lo={fill_D_bc='symmetry',
               fill_B_bc='symmetry',
               fill_psi_bc='symmetry',
               fill_phi_bc='symmetry',
               fill_ep_bc='symmetry',
               fill_mu_bc='symmetry',
              },
           hi={fill_D_bc='symmetry',
               fill_B_bc='symmetry',
               fill_psi_bc='symmetry',
               fill_phi_bc='symmetry',
               fill_ep_bc='symmetry',
               fill_mu_bc='symmetry',
              },
           },
      },
  },
}

-- === SOURCE TERMS ===
--TODO put back in first 
actions = {
  plasma = {
    type='Lorentz',
    solver   = 'implicit',
    options  = {order=3},
    states = {'ion', 'electron', 'field'}
  },

  hydro_fluxes = {
    type = 'CTU',
    corner_transport=true,
    states = {'ion', 'electron'},
  },
}

-- === GEOMETRY ===
--TODO put back in slowly once running 

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
  },
  functions = {
      charge_density=charge_density,
      div_D=div_D,
      err_div_D=err,
  },
}

