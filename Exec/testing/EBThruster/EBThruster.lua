-- ======== PROBLEM ==========
--# Assorted constants for calculations - dont change #--
--
-- e.g. 
print('Reading Lua')
mu_0_dim = 1.25663706e-6; 
ep_0_dim = 8.85418782e-12;
kb = 1.38064852e-23;
ref_lightspeed = 299792458.0; --dimensional speed of light 

--# Reference Values 
q_ref = 1.60217662e-19 -- Coulombs

m_p = 1.6726219000e-27 -- Proton mass
mfac = 1e17  --TODO problem specific for ionisation fraction or for stability? 
m_ref = mfac * m_p

n_ref = 1.0e23

-- REFERENCE VALUES 
ref_length = 0.1 -- m 
ref_mass = m_ref
ref_density = n_ref*ref_mass -- use the numver density multipled by the ref mass 


--ref_temp = 273.0
ref_temp = 1.2097998979e+06

--Larmor = 0.1
--Debye = 0.1
Larmor = 1.0e-2 --TODO google what these are/represent
Debye = 2.4e-6

phi_mag=1.1674743331935e-04 -- target D_y nondim
A_mag=1.0 -- target B_z non-dim


-- === SETTINGS ===

verbosity = 1
cfl = 0.25

do_face_sources = 0
do_CTU = 1

-- === DEFINE STATES ===

shock_x = -500

shock_mach = 2.0
density = 1.0

mass_ion = 14.0
mass_electron = 0.01

gam = 5/3
pressure = 0.5
--axis_beta = {0,0,0}

-- computations

p1 = pressure
p0 = p1*(1 + ((2*gam)/(gam+1))*(shock_mach^2 - 1))

rho1 = density
rho0 = rho1/(1 - (2/(gam+1))*(1 - 1/shock_mach^2))

a0 = math.sqrt(gam*p0/rho0)
a1 = math.sqrt(gam*p1/rho1)

u1 = 0.0
u0 = shock_mach*a1 - a0*math.sqrt(((gam-1)*shock_mach^2 + 2)/(2*gam*shock_mach^2 - (gam-1)))


-- overwrite:
p1=4.13291e-06
rho1=0.014
u1=0.03
--u1=0.0


-- functions
function shock_interface(x, L, R)
    if x <= shock_x then
	    return L
    else
	    return R
    end
end


-- fluids

function density(dat)

    x = dat['x']
    y = dat['y']

    rho = shock_interface(x, rho0, rho1)

    return rho
end

function ion_density(dat)
    return density(dat)
end

function electron_density(dat)
    return mass_electron*density(dat)/mass_ion
end

function tracer(dat)
    x = dat['x']
    y = dat['y']
    t = shock_interface(x, 0, 1)
    
    return t
end

function pressure(dat)
    x = dat['x']
    y = dat['y']
    return shock_interface(x, p0, p1)
end

function velocity_x(dat)
    x = dat['x']
    y = dat['y']
    return shock_interface(x, u0, u1)
end





states = {

  electron = {
    type='hydro',
    mass=mass_electron, 
    charge=-1.0, 
    gamma=5/3, 
    reconstruction='minmod',  
    flux='HLLC',
    refine_grad_threshold = {rho=0.1},
    value = {
      rho   = electron_density,
      x_vel = velocity_x,
      p     = pressure,
      alpha = tracer,
    },
    bc = {
            x={
                lo={
                    fill_hydro_bc = 'inflow',
                    rho=electron_density, x_vel=velocity_x, y_vel=0.0, z_vel=0.0, T=0.0041329148800799, p=pressure, alpha=tracer,
                },
            },
            },

  },

  ion = {
    type='hydro',
    mass=1.0,  
    charge= 1.0, 
    gamma=5/3, 
    reconstruction='minmod', 
    flux='HLLC',
    refine_grad_threshold = {rho=0.1},
    value = {
      rho   = ion_density,
      x_vel = velocity_x,
      p     = pressure,
      alpha = tracer,
    },
    bc = {
            x={
                lo={
                    fill_hydro_bc = 'inflow',
                    rho=ion_density, x_vel=velocity_x, y_vel=0.0, z_vel=0.0, T=0.0041329148800799, p=pressure, alpha=tracer,
                },
            },
            },


  },

  field = {
    type='field',
    static=1,
    reconstruction='O6',
    value={x_D=0.0,y_D=phi_mag,z_D=0.0,x_B=0.0,y_B=0.0,z_B=A_mag},
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

-- === SOURCES ===

sources = {

  plasma={
      solver = 'implicit',
      options = {order=3},
      sources = {
          plasma={'ion', 'electron', 'field',
              type='Lorentz',
          },
      },
  },    
}


-- === GEOMETRY ===

-- options
refine_cutcells = 1
merge_fraction = 0.5

wire_length = 3
wire_thickness = 0.2
wire_separation = 1.6
shield_radius = 0.25

Az_magnitude = 1.0

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

-- collection = {{{x_lo, x_hi},{y_lo,y_hi}}, ... }
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

function shield_def(x,y)
  local circles = {
    {shield_radius, 0.0, 0.0},
  }
  return make_circles(x,y,circles)
end


function wire_1(x,y)
  local rect = {
    {{-wire_length/2, wire_length/2},{wire_separation, wire_separation+wire_thickness}},
  }
  return make_rectangles(x,y,rect)
end

function wire_2(x,y)
  local rect = {
    {{-wire_length/2, wire_length/2},{-wire_separation-wire_thickness, -wire_separation}},
  }
  return make_rectangles(x,y,rect)
end

function topblock(x,y)
  local rect = {
    {{-5, 5},{1.5, 3}},
  }
  return make_rectangles(x,y,rect)
end

function botblock(x,y)
  local rect = {
    {{-5, 5},{-3, -1.5}},
  }
  return make_rectangles(x,y,rect)
end

function outerbound(x,y)
  local rect = {
    {{-4, -2.9},{-3, 3}},
  }
  return make_rectangles(x,y,rect)
end

function outerboundhi(x,y)
  local rect = {
    {{2.9,4},{-3, 3}},
  }
  return make_rectangles(x,y,rect)
end


function phi_f(dat)
  local phival
  x = dat['x']
  y = dat['y']

  phival=y*phi_mag
  return phival
end

function phi_fm(dat)
  local phival
  x = dat['x']
  y = dat['y']

  phival=-y*phi_mag
  return phival
end


function Ay_f(dat)
  local Aval
  x = dat['x']
  y = dat['y']

  Aval=x*A_mag
  return Aval
end

embedded_boundaries = {

  --[[shield = {
    geom=shield_def,
    bcs={
      ion={
        type='slip_wall'
      },
      electron={
        type='slip_wall'
      },
    },
    boolean_operation='and',
    inside=0,
  },--]]

  topwall = {
    geom=topblock,
    bcs={
      ion={
        type='slip_wall'
      },
      electron={
        type='slip_wall'
      },
    },
    boolean_operation='and',
    inside=0,
  },

  botwall = {
    geom=botblock,
    bcs={
      ion={
        type='slip_wall'
      },
      electron={
        type='slip_wall'
      },
    },
    boolean_operation='and',
    inside=0,
  },


  solenoid_part_1 = {
    geom=wire_1,
    bcs={
      field={
        types={
          'scalar_potential', 'vector_potential'
        }, 
        phi=-phi_mag*wire_separation, 
        A0=0.0,
        A1=Ay_f,
        A2=0.0,
        align_with_boundary=false,
      },
    },
    boolean_operation='or',
    inside=0,
  },

  solenoid_part_2 = {
    geom=wire_2,
    bcs={
      field={
        types={
          'scalar_potential', 'vector_potential'
        }, 
        phi=phi_mag*wire_separation, 
        A0=0.0,
        A1=Ay_f,
        A2=0.0,
        align_with_boundary=false,
      },
    },
    boolean_operation='and',
    inside=0,
  },
  
    --[[zBgenlo = {
    geom=outerbound,
    bcs={
      field={
        types={
          'scalar_potential', 'vector_potential'
        }, 
        phi=phi_f, 
        A0=0.0,
        A1=Ay_f,
        A2=0.0,
        align_with_boundary=false,
      },
    },
    boolean_operation='and',
    inside=0,
  },

    zBgenhi = {
    geom=outerboundhi,
    bcs={
      field={
        types={
          'scalar_potential', 'vector_potential'
        }, 
        phi=phi_fm, 
        A0=0.0,
        A1=Ay_f,
        A2=0.0,
        align_with_boundary=false,
      },
    },
    boolean_operation='and',
    inside=0,
  },--]]

  
  --[[solenoid_part_1 = {
    geom=wire_1,
    bcs={
      field={
        types={
          'surface_charge', 'surface_current'
        }, 
        charge=1.0, 
        j1=0.0,
       },
    },
    boolean_operation='or',
    inside=0,
  },

  solenoid_part_2 = {
    geom=wire_2,
    bcs={
      field={
        types={
          'surface_charge', 'surface_current'
        }, 
        charge=-1.0, 
        j1=0.0,
       },
    },
    boolean_operation='and',
    inside=0,
  }--]]
  
  
}

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
print("Read Lua")

