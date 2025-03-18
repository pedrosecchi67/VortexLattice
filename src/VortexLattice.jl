module VortexLattice

  using LinearAlgebra
  using DelimitedFiles

  include("MeshPlotting.jl")
  using .MeshPlotting

  """
  Geometric tolerance
  """
  geps=1e-18

  """
  Function to define inf. velocity of a horseshoe vortex

  * `a`: first kin (relative to collocation point)
  * `b`: second kin (relative to collocation point)
  * `x`: wind direction (unit vector)
  
  * return: the influence velocity over the collocation point
  """
  function hshoe_vinf(a::Vector{Float64}, b::Vector{Float64}, x::Vector{Float64})
    na=norm(a)
    nb=norm(b)
    ab=dot(a, b)

    d1=(na*nb+ab)
    if d1>-geps && d1<=0.0
      d1=-geps
    elseif d1<geps && d1>0.0
      d1=geps
    end

    d2=(na-dot(a, x))
    if d2>-geps && d2<=0.0
      d2=-geps
    elseif d2<geps && d2>0.0
      d2=geps
    end

    d3=(nb-dot(b, x))
    if d3>-geps && d3<=0.0
      d3=-geps
    elseif d3<geps && d3>0.0
      d3=geps
    end

    vbound = (cross(a, b)/d1)*(1.0/na+1.0/nb)
    vinf = cross(a, x)/(na*d2).-cross(b, x)/(nb*d3)

    vinf./=(4*pi)
    vbound./=(4*pi)

    return vinf, vbound
  end

  """
  Function to define inf. velocity of a finite vortex line segment

  * `a`: first kin (relative to collocation point)
  * `b`: second kin (relative to collocation point)
  
  * return: the influence velocity over the collocation point
  """
  function vline_segment_vinf(a::Vector{Float64}, b::Vector{Float64})
    na=norm(a)
    nb=norm(b)
    ab=dot(a, b)

    d1=(na*nb+ab)
    if d1>-geps && d1<=0.0
      d1=-geps
    elseif d1<geps && d1>0.0
      d1=geps
    end

    vinf=(cross(a, b)/d1)*(1.0/na+1.0/nb)

    vinf/=(4*pi)

    return vinf
  end

  """
  Struct to contain information about an airfoil

  * `camberline_coords`: matrix `(shape (npts, 2))` containing airfoil camber line coordinates
  * `CL_correct`: function as ``C_L=C_{L, inv}`` `+CL_correct(CL)`
  * `CD_correct`: function as ``C_D=C_{D, inv}``+ `CD_correct(CL)`
  """
  mutable struct Airfoil
    camberline_coords::Matrix{Float64}
    CL_correct::Function
    CD_correct::Function
  end

  """
  Constructor for inviscid flat plate airfoil
  """
  inviscid_flatplate()=Airfoil(reshape([0.0, 1.0, 0.0, 0.0], 2, 2), CL->0.0, CL->0.0)

  function _interp_intra(afl_pts::Matrix{Float64}, x::Float64)

    for i = size(afl_pts, 1):-1:2
      if afl_pts[i, 1] <= x && afl_pts[i - 1, 1] >= x
        return afl_pts[i, 2] + (afl_pts[i - 1, 2] - afl_pts[i, 2]) * (x - afl_pts[i, 1]) / (afl_pts[i - 1, 1] - afl_pts[i, 1])
      end
    end

    throw(
      error(
        "Point at $x not found in airfoil"
      )
    )

  end

  """
  Constructor for airfoil based on AVL correction coefficients

  * `airfoil_file`: file for airfol in Selig format
  * `CLAF`: ratio between lift coefficient derivative `CLα` and `2π`, for lift correction
  * `polar_points`: vector of tuples with points for polar interpolation.

  Example:

  ```
  [
    (0.0, 0.012),
    (0.2, 0.010),
    (0.4, 0.012) # arbitraty parabolic interpolation from `CL` to `CD`
  ]
  ```
  """
  function Airfoil(
    afl_file::String,
    CLAF::Float64,
    polar_points::Vector{Tuple{Float64, Float64}}
  )

    try
      airfoil_pts = readdlm(afl_file; skipstart = 1)

      nLE = argmin(airfoil_pts[:, 1])

      upper_surface_points = airfoil_pts[nLE:-1:1, :]

      for i = 1:nLE
        upper_surface_points[i, 2] = (upper_surface_points[i, 2] + _interp_intra(airfoil_pts, upper_surface_points[i, 1])) / 2
      end

      CLs = [
        CL for (CL, _) in polar_points
      ]
      CDs = [
        CD for (_, CD) in polar_points
      ]

      pA = [
        CLs.^0 CLs.^1 CLs.^2
      ]

      pcoefs = pA \ CDs

      return Airfoil(
        upper_surface_points,
        CL -> (CLAF - 1.0) * CL,
        CL -> (pcoefs[1] + pcoefs[2] * CL + pcoefs[3] * CL ^ 2)
      )
    catch
      throw(
        error(
          "Error reading from file $afl_file"
        )
      )
    end

  end

  """
  Function to get angle of attack due to ``dz/dx`` component of airfoil geometry
  
  * `afl`: the airfoil
  * `bdisc_colpts`: stations identifying positions of collocation points as a percentage of the chord
  
  * return: angles at each chordwise station
  """
  function get_alpha0s(afl::Airfoil, bdisc_colpts::Vector{Float64})
    npts=size(bdisc_colpts, 1)
    nafl=size(afl.camberline_coords, 1)

    alphas=zeros(Float64, npts)

    for i=1:npts
      for j=1:(nafl-1)
        if afl.camberline_coords[j, 1]<bdisc_colpts[i, 1] && afl.camberline_coords[j+1, 1]>bdisc_colpts[i, 1]
          alphas[i]=atan((afl.camberline_coords[j+1, 2]-afl.camberline_coords[j, 2]), (afl.camberline_coords[j+1, 1]-afl.camberline_coords[j, 1]))

          break
        end
      end
    end

    return alphas
  end

  export Airfoil, inviscid_flatplate

  """
  Struct to define wing sections

  * `LE`: position of the leading edge
  * `chord`: section chord
  * `incidence`: section incidences
  * `afl`: instance of struct Airfoil to use for viscous corrections
  """
  mutable struct Section
    LE::Vector{Float64}
    chord::Float64
    incidence::Float64
    afl::Airfoil
    controls::Vector{Tuple{Symbol, Float64, Float64}}
  end

  """
  Constructor for wing section. Uses a flat plate airfoil as default

  * `LE`: position of the leading edge
  * `chord`: section chord
  * `incidence`: section incidence. Defaults to zero
  * `afl`: instance of struct Airfoil to use for viscous corrections. Defaults to flat plate
  * `controls`: array of tuples indicating present control surfaces. Tuples should indicate, respectively, the control's name (Symbol), the percentage of the chord in
    which it's initiated and its gain in the present section
  """
  Section(LE::Vector{Float64}, chord::Float64; incidence::Float64=0.0, afl=inviscid_flatplate(), controls=[])=Section(LE, chord, incidence, afl, controls)

  export Section

  """
  Struct to obtain a quadrant between two wing sections

  * `sec_left`: left section
  * `sec_right`: right section
  * `mtosys`: matrix for local coordinate system (keeps x axis unchanged, y and z are altered so that y is parallel to the surface)
  """
  mutable struct Quadrant
    sec_left::Section
    sec_right::Section
    mtosys::Matrix{Float64}
  end

  export Quadrant

  """
  Struct to model a lifting surface

  * `quadrants`: vector of wing quadrants
  * `pts`: points representing panel corners `(shape (nm+1, nn+1, 3))`
  * `colpts`: locations of collocation points `(shape (nm, nn, 3))`
  * `normals`: normal vectors at each collocation point `(shape (nm, nn, 3))`
  * `mtosys`: local coordinate system matrix at each panel
  * `kins`: locations for vortex kins `(shape (nm, nn+1, 3))`
  * `chords`: chord at each panel strip
  * `widths`: width of each panel strip
  * `quarter_chords`: matrix `(shape (nstrips, 3))` with quarter chord positions across the span
  * `airfoil_ponderations`: array of tuples containing two Airfoil instances and a ponderation factor. Used for viscous corrections
  * `control_ponderations`: matrix of tuples indicating, respectively, the name (Symbol) and the gain (float) for each control surface acting on each panel
  * `mac`: mean aerodynamic chord
  * `S`: surface area
  """
  mutable struct Surface
    quadrants::Vector{Quadrant}
    pts::Array{Float64, 3}
    colpts::Array{Float64, 3}
    normals::Array{Float64, 3}
    mtosys::Matrix{Matrix{Float64}}
    kins::Array{Float64, 3}
    inds::Matrix{Int64}
    chords::Vector{Float64}
    widths::Vector{Float64}
    quarter_chords::Matrix{Float64}
    airfoil_ponderations::Vector{Tuple{Float64, Airfoil, Airfoil}}
    control_ponderations::Matrix{Vector{Tuple{Symbol, Float64, Float64}}}
    mac::Float64
    S::Float64
  end

  """
  Function returning whether or not a section possesses a similar control surface definition

  * `sect`: the section instance in which to check
  * `cont`: tuple with, respectively, the control surface's name (Symbol), its chord percentage in the section and its gain as present in Section.controls

  * return: an index for the similar control surface in sect, if present; 0 if otherwise
  """
  function has_similar_control(sect::Section, cont::Tuple{Symbol, Float64, Float64})
    namesymb, _, gain=cont

    for (i, compcont) in enumerate(sect.controls)
      if namesymb==compcont[1] && gain==compcont[3]
        return i
      end
    end

    return 0
  end

  """
  Function to obtain a cossenoidal distribution between 0 and 1, with `n` points
  """
  function cossenoidal_distribution(n::Int64)

    (sin.(collect(LinRange(- pi / 2, pi / 2, n))) .+ 1.0) ./ 2

  end

  export cossenoidal_distribution

  """
  Constructor for a lifting surface based on a list of sections (from left to right)

  For arguments `cdiscs` and `bdiscs`, a usage example is:

  `sects=[
    (LE=[1.0, -1.0, 0.0], chord=0.5),
    (LE=[0.0, 0.0, 0.0], chord=1.0),
    (LE=[1.0, -1.0, 0.0], chord=0.5)
  ]
  bdiscs=[
    [0.0, 0.5, 1.0],
    [0.0, 1.0]
  ]
  cdiscs=[
    [0.0, 0.5, 1.0],
    [0.0, 0.5, 1.0],
    [0.0, 0.5, 1.0]
  ]

  would result in the panel corners:
  [
    [1.0, -1.0, 0.0] [0.5, -0.5, 0.0] [0.0, 0.0, 0.0] [1.0, 1.0, 0.0]
    [1.25, -1.0, 0.0] [0.875, -0.5, 0.0] [0.5, 0.0, 0.0] [1.25, 1.0, 0.0]
    [1.5, -1.0, 0.0] [1.25, -0.5, 0.0] [1.0, 0.0, 0.0] [1.5, 1.0, 0.0]
  ]`

  * `sects`: vector of sections to base the surface on
  * `cdiscs`: vector of vectors, each describing a mapping for the chordwise discretization in a wing section (ranging from zero to 1, and spaced as desired for panel
    edges)
  * `bdiscs`: vector of vectors, each describing a mapping for the spanwise discretization in a wing quadrant
  """
  function Surface(sects::Vector{Section}, cdiscs::Vector{Vector{Float64}}, bdiscs::Vector{Vector{Float64}})
    ns=size(sects, 1)
    
    if size(cdiscs, 1)!=ns
      throw(error("Surface:ERROR:Inconsistent number of chordwise discretization definitions provided for sections ($ns sections expected)"))
    end
    
    if size(bdiscs, 1)!=ns-1
      throw(error("Surface:ERROR:nconsistent number of spanwise discretization definitions provided for wing quadrants ($(ns-1) sections expected)"))
    end

    quads=Vector{Quadrant}()

    for i=1:(size(sects, 1)-1)
      toright=sects[i+1].LE.-sects[i].LE
      toright[1]=0.0
      toright/=norm(toright)

      mtosys=zeros(Float64, 3, 3)

      mtosys[1, 1]=1.0
      mtosys[2, :]=toright
      mtosys[3, 2]=-toright[3]
      mtosys[3, 3]=toright[2]

      push!(quads, Quadrant(sects[i], sects[i+1], mtosys))
    end

    nbs=sum([size(bd, 1)-1 for bd in bdiscs])+1
    ncs=size(cdiscs[1], 1)

    for i=2:size(cdiscs, 1)
      if size(cdiscs[i], 1)!=ncs
        throw(error("Surface:ERROR:wing section $i has been indicated a different number of chordwise panels when put side by side with section 1"))
      end
    end

    pts=zeros(Float64, ncs, nbs, 3)
    normals=zeros(Float64, ncs-1, nbs-1, 3)
    panel_mtosys=Matrix{Matrix{Float64}}(undef, ncs-1, nbs-1)

    airfoil_ponderations=Vector{Tuple{Float64, Airfoil, Airfoil}}()

    control_ponderations=Matrix{Vector{Tuple{Symbol, Float64, Float64}}}(undef, ncs-1, nbs-1)
    for i=1:ncs-1
      for j=1:nbs-1
        control_ponderations[i, j]=Vector{Tuple{Symbol, Float64, Float64}}()
      end
    end

    nb=1

    for i=1:(ns-1)
      for j=1:(i==ns-1 ? size(bdiscs[i], 1) : size(bdiscs[i], 1)-1)
        factor=bdiscs[i][j]

        cdisc_med=factor*cdiscs[i+1][end:-1:1].+(1.0-factor)*cdiscs[i][end:-1:1]

        chord_med=factor*sects[i+1].chord+(1.0-factor)*sects[i].chord
        inc_med=deg2rad(factor*sects[i+1].incidence+(1.0-factor)*sects[i].incidence)

        pts[:, nb, 1].=(chord_med*cdisc_med).+(factor*sects[i+1].LE[1]+(1.0-factor)*sects[i].LE[1])
        pts[:, nb, 2].=factor*sects[i+1].LE[2]+(1.0-factor)*sects[i].LE[2]
        pts[:, nb, 3].=factor*sects[i+1].LE[3]+(1.0-factor)*sects[i].LE[3]

        if nb<nbs
          cdisc_med=(cdisc_med[1:(end-1)].+cdisc_med[2:end]*3)/4
          alphas_med=factor*get_alpha0s(sects[i+1].afl, cdisc_med).+get_alpha0s(sects[i].afl, cdisc_med)*(1.0-factor)
          sins=[sin(a-inc_med) for a in alphas_med]
          coss=[cos(a-inc_med) for a in alphas_med]

          normals[:, nb, 1].=-sins
          normals[:, nb, 2].=coss.*quads[i].mtosys[3, 2]
          normals[:, nb, 3].=coss.*quads[i].mtosys[3, 3]

          for j=1:ncs-1
            panel_mtosys[j, nb]=quads[i].mtosys
          end

          push!(airfoil_ponderations, (factor, sects[i].afl, sects[i+1].afl))

          for cont in sects[i].controls
            sim=has_similar_control(sects[i+1], cont)

            if sim!=0
              symb, perc1, gain=cont
              _, perc2, _=sects[i+1].controls[sim]

              perc=1.0-(factor*perc2+(1.0-factor)*perc1)

              for j=1:ncs-1
                arm=cdisc_med[j]-perc

                if arm>0.0
                  push!(control_ponderations[j, nb], (symb, arm, gain))
                end
              end
            end
          end
        end

        nb+=1
      end
    end

    chords=zeros(Float64, nbs-1)
    widths=zeros(Float64, nbs-1)

    for i=1:(nbs-1)
      chords[i]=((pts[1, i, 1]-pts[end, i, 1])+(pts[1, i+1, 1]-pts[end, i+1, 1]))/2
      widths[i]=sqrt((pts[1, i, 2]-pts[1, i+1, 2])^2+(pts[1, i, 3]-pts[1, i+1, 3])^2)

      for j=1:(ncs-1)
        for k=1:size(control_ponderations[j, i], 1)
          symb, arm, gain=control_ponderations[j, i][k]
          arm*=chords[i]

          control_ponderations[j, i][k]=(symb, arm, gain)
        end
      end
    end

    quarter_chords=((pts[end, 1:(end-1), :]*1.5).+(pts[end, 2:end, :]*1.5).+(pts[1, 1:(end-1), :]*0.5).+(pts[1, 2:end, :]*0.5))/4

    colpts=(1.5*pts[1:(end-1), 1:(end-1), :].+1.5*pts[1:(end-1), 2:end, :].+0.5*pts[2:end, 1:(end-1), :].+0.5*pts[2:end, 2:end, :])/4
    kins=(pts[1:(end-1), :, :].+3.0*pts[2:end, :, :])/4

    S=sum(chords.*widths)
    mac=sum((chords.^2).*widths)/S

    return Surface(quads, pts, colpts, normals, panel_mtosys, kins, zeros(Int64, ncs-1, nbs-1), 
        chords, widths, quarter_chords, airfoil_ponderations, control_ponderations, 
        mac, S)
  end

  export Surface

  """
  Abstract data type to contain information about an aircraft

  * `Sref`: reference surface
  * `cref`: reference chord
  * `bref`: reference span
  * `surfaces`: array of surfaces to compose the aircraft
  * `colpts`: the totality of the mesh's control points (vector of vectors format)
  * `normals`: their respective normal vectors (vector of vectors format)
  * `mtosys`: vector of local coordinate system matrices
  * `control_ponderations`: vector of vector of tuples indicating, respectively, the name (Symbol), 
    the hinge arm and the gain for each control surface acting on each panel
  * `kins`: their respective vortex kins (vector of tuple of vectors format)
  * `AICM3`: aerodynamic influence coefficient matrix in 3D
  * `downwash_AICM3`: aerodynamic influence coefficient matrix in 3D disregarding bound vortex segments 
    (check Katz & Plotkin for more information)
  * `AICM`: aerodynamic influence coefficient matrix
  * `LU_AICM`: LU-decomposed aerodynamic coefficient matrix
  """
  mutable struct Aircraft
    Sref::Float64
    cref::Float64
    bref::Float64
    surfaces::Vector{Surface}
    colpts::Vector{Vector{Float64}}
    normals::Vector{Vector{Float64}}
    mtosys::Vector{Matrix{Float64}}
    control_ponderations::Vector{Vector{Tuple{Symbol, Float64, Float64}}}
    kins::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
    AICM3::Array{Float64, 3}
    downwash_AICM3::Array{Float64, 3}
    AICM::Matrix{Float64}
    LU_AICM::LU
  end

  """
  Constructor for Aircraft struct based on an array of surfaces

  * `surfaces`: vector of surfaces
  * `Sref`: reference surface
  * `cref`: reference chord
  * `bref`: reference span
  * `ground_effect`: boolean for ground effect analysis
  * `ground_height`: ground distance from main surface
  """
  function Aircraft(surfaces::Vector{Surface}; Sref::Float64=1.0, cref::Float64=1.0, bref::Float64=1.0,ground_effect::Bool=false,ground_height::Float64=0.0)
    colpts=Vector{Vector{Float64}}()
    normals=Vector{Vector{Float64}}()
    mtosys=Vector{Matrix{Float64}}()
    kins=Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
    control_ponderations=Vector{Vector{Tuple{Symbol, Float64, Float64}}}()

    ncol=0
    for surf in surfaces
      nm, nn, _=size(surf.colpts)

      for i=1:nm
        for j=1:nn
          ncol+=1

          push!(colpts, surf.colpts[i, j, :])
          push!(normals, surf.normals[i, j, :])
          push!(mtosys, surf.mtosys[i, j])
          push!(control_ponderations, surf.control_ponderations[i, j])
          push!(kins, (surf.kins[i, j, :], surf.kins[i, j+1, :]))

          surf.inds[i, j]=ncol
        end
      end
    end

    AICM3=zeros(Float64, ncol, ncol, 3)
    downwash_AICM3=zeros(Float64, ncol, ncol, 3)
    AICM=zeros(Float64, ncol, ncol)

    xhat=[1.0, 0.0, 0.0]

    for (i, (colpt, normal)) in enumerate(zip(colpts, normals))
      for (j, kin) in enumerate(kins)
        a=colpt.-kin[1]
        b=colpt.-kin[2]
        vinf, vbound=hshoe_vinf(a, b, xhat)
        vinfluence = vinf .+ vbound

        # ground effect condition
        if ground_effect
          # reflecting points to create mirror image
          kin1_m = [kin[1][1], kin[1][2], -2*ground_height - kin[1][3]]
          kin2_m = [kin[2][1], kin[2][2], -2*ground_height - kin[2][3]]
          a_m = colpt .- kin1_m
          b_m = colpt .- kin2_m
          vinf_m, vbound_m = hshoe_vinf(a_m, b_m, xhat)
          # mirror vortice influence
          mirror_influence = -(vinf_m .+ vbound_m)
          vinfluence .+= mirror_influence
        end

        AICM3[i, j, :].=vinfluence
        downwash_AICM3[i, j, :].=vinf
        AICM[i, j]=dot(vinfluence, normal)
      end
    end

    LU_AICM=lu(AICM)
    
    return Aircraft(Sref, cref, bref, surfaces, colpts, normals, mtosys, control_ponderations, kins, AICM3, downwash_AICM3, AICM, LU_AICM)
  end

  export Aircraft

  """
  Function to obtain air velocities at each mesh collocation point

  * `acft`: the aircraft object
  * `alpha`: angle of attack (degrees)
  * `beta`: sideslip angle (degrees)
  * `p`: x axis (stability and control coordinate system, rad) angular velocity
  * `q`: y axis (stability and control coordinate system, rad) angular velocity
  * `r`: z axis (stability and control coordinate system, rad) angular velocity

  * return: a matrix `(shape (nvortexes, 3))` with the velocity at each collocation point
  """
  function get_freestream_velocities(acft::Aircraft; alpha::Float64=0.0, beta::Float64=0.0, p::Float64=0.0, q::Float64=0.0, r::Float64=0.0)
    omega=[
      p*2.0/acft.bref,
      q*2.0/acft.cref,
      r*2.0/acft.bref
    ]

    xhat=[cosd(alpha)*cosd(beta), cosd(alpha)*sind(beta), sind(alpha)]

    npts=size(acft.colpts, 1)

    vels=zeros(Float64, npts, 3)

    for (i, colpt) in enumerate(acft.colpts)
      vels[i, :].=xhat.-cross(omega, colpt)
    end

    return vels
  end

  """
  Function to obtain rotation matrix for a given control surface's deflection

  * `angle`: angle of deflection
  * `mtosys`: matrix for conversion to local coordinate system

  * return: a matrix to multiply the normal vector by so as to correctly rotate it
  """
  function get_rotation_matrix(angle::Float64, mtosys::Matrix{Float64})
    rot=Matrix{Float64}(undef, 3, 3)

    rot[1, 1]=cosd(angle)
    rot[1, 2]=0.0
    rot[1, 3]=sind(angle)

    rot[2, 1]=0.0
    rot[2, 2]=1.0
    rot[2, 3]=0.0

    rot[3, 1]=-sind(angle)
    rot[3, 2]=0.0
    rot[3, 3]=cosd(angle)

    return transpose(mtosys)*rot*mtosys
  end

  """
  Function to return an array of vortex intensities comprising a vortex lattice method solution for the given flight conditions
  
  * `acft`: the aircraft object
  * `alpha`: angle of attack
  * `beta`: sideslip angle
  * `p`: x axis (stability and control coordinate system) angular velocity
  * `q`: y axis (stability and control coordinate system) angular velocity
  * `r`: z axis (stability and control coordinate system) angular velocity
  * `control_deflections`: dictionary relating control naming (Symbols) to control deflections

  * return: an array with the vorticity strength in each horseshoe vortex in the mesh
  * return: total (influence+freestream) velocity at each collocation point
  """
  function VLsolve(acft::Aircraft; alpha::Float64=0.0, beta::Float64=0.0, p::Float64=0.0, q::Float64=0.0, r::Float64=0.0, 
      control_deflections::Dict{Symbol, Float64}=Dict{Symbol, Float64}())
    vels=get_freestream_velocities(acft; alpha=alpha, beta=beta, p=p, q=q, r=r)

    npts=size(acft.colpts, 1)

    normal_vels=Vector{Float64}(undef, npts)

    for i=1:npts
      conts=acft.control_ponderations[i]

      angle=0.0
      for (contname, arm, gain) in conts
        if haskey(control_deflections, contname)
          angle+=control_deflections[contname]*gain
        end
      end

      normal_vels[i]=dot(get_rotation_matrix(angle, acft.mtosys[i])*acft.normals[i], vels[i, :])
    end

    soln=-(acft.LU_AICM\normal_vels)

    for i=1:3
      vels[:, i].+=acft.downwash_AICM3[:, :, i]*soln
    end

    return soln, vels
  end

  export VLsolve

  """
  Plot aircraft geometry using matplotlib backend
  
  * `acft`: the aircraft to plot
  """
  plot_aircraft(acft::Aircraft)=plot_grids([s.pts for s in acft.surfaces])

  export plot_aircraft

  """
  Plot a property throughout an aircraft's mesh

  * `acft`: the aircraft at hand
  * `prop`: array mapping the property's value at every panel
  """
  plot_props(acft::Aircraft, prop::Vector{Float64})=plot_props_in_mesh(
      [s.pts for s in acft.surfaces], [prop[s.inds] for s in acft.surfaces])

  export plot_props

  """
  Obtain forces on each panel based on past solution

  * `soln`: vorticity vector as returned from VLsolve
  * `vels`: velocity vector retured from VLsolve

  * return: an array with forces at each panel, with the same shape as vels
  """
  function get_forces(acft::Aircraft, soln::Vector{Float64}, vels::Matrix{Float64})
    forces=zeros(Float64, size(soln, 1), 3)

    for (i, (vort, kin)) in enumerate(zip(soln, acft.kins))
      forces[i, :].=cross(vels[i, :], kin[2].-kin[1])*(2*vort)
    end

    return forces
  end

  export get_forces

  """
  Get arrays with forces and moments around the quarter chord at each panel strip

  * `forces`: array of inviscid forces, as obtained from get_forces
  * `acft`: Aircraft instance
  * `surf`: Surface instance
  * `alpha`: angle of attack
  * `beta`: sideslip angle

  * return: a dictionary with coefficient arrays (keys: ``:CLi; :CL; :CDi; :CD; :CX; :CY; :CZ; :Cl; :Cm; :Cn`). Moments are obtained around the quarter chord,
      And refer to the local chord (check variables surf.quarter_chords and surf.chords)
  """
  function get_strip_coefficients(forces::Matrix{Float64}, surf::Surface; alpha::Float64=0.0, beta::Float64=0.0)
    nstrips=size(surf.chords, 1)
    nc=size(surf.inds, 1)

    xhat=[cosd(alpha)*cosd(beta), cosd(alpha)*sind(beta), sind(alpha)]
    zhat=[-sind(alpha)*cosd(beta), -sind(alpha)*sind(beta), cosd(alpha)]

    strip_forces=Dict(
        :CLi=>zeros(nstrips),
        :CL=>zeros(nstrips),
        :CDi=>zeros(nstrips),
        :CD=>zeros(nstrips),
        :CX=>zeros(nstrips),
        :CY=>zeros(nstrips),
        :CZ=>zeros(nstrips),
        :Cl=>zeros(nstrips),
        :Cm=>zeros(nstrips),
        :Cn=>zeros(nstrips)
      )

    for i=1:nstrips
      force=zeros(Float64, 3)
      mom=zeros(Float64, 3)

      for j=1:nc
        locforc=forces[surf.inds[j, i], :]

        force.+=locforc
        mom.+=cross(surf.colpts[j, i, :].-surf.quarter_chords[i, :], locforc)
      end

      chord=surf.chords[i]
      area=chord*surf.widths[i]

      force/=area
      mom/=(area*chord)

      strip_forces[:CX][i]=force[1]
      strip_forces[:CY][i]=force[2]
      strip_forces[:CZ][i]=force[3]

      strip_forces[:Cl][i]=mom[1]
      strip_forces[:Cm][i]=mom[2]
      strip_forces[:Cn][i]=mom[3]

      CDi=dot(xhat, force)
      CLi=dot(zhat, force)

      strip_forces[:CDi][i]=CDi
      strip_forces[:CLi][i]=CLi

      factor, afl1, afl2=surf.airfoil_ponderations[i]

      CL=CLi+(1.0-factor)*afl1.CL_correct(CLi)+factor*afl2.CL_correct(CLi)
      CD=CDi+(1.0-factor)*afl1.CD_correct(CLi)+factor*afl2.CD_correct(CLi)

      strip_forces[:CL][i]=CL
      strip_forces[:CD][i]=CD
    end

    return strip_forces
  end

  export get_strip_coefficients

  """
  Function to obtain individual surface forces

  * `surf`: Surface instance
  * `forces`: forces as output from get_forces
  * `vels`: velocities as output from get_freestream_velocities
  * `alpha`: angle of attack
  * `beta`: angle of sideslip
  * `CG`: center of gravity around which to calculate moments

  * return: strip forces (dictionary with keys) as output from `get_strip_coefficients`
  * return: surface forces: a dictionary with the same keys, but considering total forces on the surface, with reference to its own surface area
    and mean aerodynamic chord (see variables `surf.S` and `surf.mac`, respectively), and moments centered around the given center of gravity
  """
  function get_surf_coefficients(surf::Surface, forces::Matrix{Float64}, vels::Matrix{Float64}; 
      alpha::Float64=0.0, beta::Float64=0.0, CG::Vector{Float64}=[0.0, 0.0, 0.0])
    strip_forces=get_strip_coefficients(forces, surf; alpha=alpha, beta=beta)

    nstrips=size(surf.chords, 1)

    xhat=[cosd(alpha)*cosd(beta), cosd(alpha)*sind(beta), sind(alpha)]
    zhat=[-sind(alpha)*cosd(beta), -sind(alpha)*sind(beta), cosd(alpha)]

    forc=zeros(Float64, 3)
    mom=zeros(Float64, 3)

    inviscid_forc=zeros(Float64, 3)
    inviscid_mom=zeros(Float64, 3)

    locforc=zeros(Float64, 3)
    locmom=zeros(Float64, 3)

    loc_inviscid_forc=zeros(Float64, 3)
    loc_inviscid_mom=zeros(Float64, 3)

    for i=1:nstrips
      cloc=surf.chords[i]
      dS=surf.widths[i]*cloc

      application_point=surf.quarter_chords[i, :]

      arm=application_point.-CG

      loc_inviscid_forc[1]=strip_forces[:CX][i]*dS
      loc_inviscid_forc[2]=strip_forces[:CY][i]*dS
      loc_inviscid_forc[3]=strip_forces[:CZ][i]*dS

      loc_inviscid_mom[1]=strip_forces[:Cl][i]*dS*cloc
      loc_inviscid_mom[2]=strip_forces[:Cm][i]*dS*cloc
      loc_inviscid_mom[3]=strip_forces[:Cn][i]*dS*cloc

      loc_inviscid_mom.+=cross(arm, loc_inviscid_forc)
      
      locforc.=loc_inviscid_forc.+(xhat*(strip_forces[:CD][i]-strip_forces[:CDi][i]).+zhat*(strip_forces[:CL][i]-strip_forces[:CLi][i]))*dS
      locmom.=loc_inviscid_mom.+cross(arm, locforc.-loc_inviscid_forc)

      forc.+=locforc
      mom.+=locmom

      inviscid_forc.+=loc_inviscid_forc
      inviscid_mom.+=loc_inviscid_mom
    end

    CL=dot(zhat, forc)/surf.S
    CD=dot(xhat, forc)/surf.S

    CLi=dot(zhat, inviscid_forc)/surf.S
    CDi=dot(xhat, inviscid_forc)/surf.S

    CX=forc[1]/surf.S
    CY=forc[2]/surf.S
    CZ=forc[3]/surf.S

    Cl=mom[1]/(surf.S*surf.mac)
    Cm=mom[2]/(surf.S*surf.mac)
    Cn=mom[3]/(surf.S*surf.mac)

    surf_forces=Dict(
          :CL=>CL,
          :CD=>CD,
          :CLi=>CLi,
          :CDi=>CDi,
          :CX=>CX,
          :CY=>CY,
          :CZ=>CZ,
          :Cl=>Cl,
          :Cm=>Cm,
          :Cn=>Cn
        )

    return strip_forces, surf_forces
  end

  """
  Function to obtain hinge moment coefficients (in reference to the aircraft's reference surface and chord).
  Hinge moments are multiplied by control gains where any is present

  * `acft`: Aircraft instance
  * `forces`: forces as output from get_forces

  * return: a dictionary with the hinge moment coefficient for each symbol present across sections
  """
  function get_hinge_moments(acft::Aircraft, forces::Matrix{Float64})
    hinge_moments=Dict{Symbol, Float64}()

    for i=1:size(acft.colpts, 1)
      for cont in acft.control_ponderations[i]
        symb, arm, gain=cont

        if !haskey(hinge_moments, symb)
          hinge_moments[symb]=0.0
        end

        hinge_moments[symb]+=arm*dot(acft.normals[i], forces[i, :])*gain
      end
    end

    for (k, v) in hinge_moments
      hinge_moments[k]=v/(acft.Sref*acft.cref)
    end

    return hinge_moments
  end

  """
  Function to obtain all pertinent coefficients for an aircraft

  * `acft`: `Aircraft` struct instance
  * `alpha`: angle of attack
  * `beta`: sideslip angle
  * `p`: x axis (stability and control coordinate system) angular velocity
  * `q`: y axis (stability and control coordinate system) angular velocity
  * `r`: z axis (stability and control coordinate system) angular velocity
  * `control_deflections`: dictionary relating control symbols and control deflections
  * `CG`: vector with coordinates for momentum calculations
  """
  function get_data(acft::Aircraft; alpha::Float64=0.0, beta::Float64=0.0,
      p::Float64=0.0, q::Float64=0.0, r::Float64=0.0, control_deflections::Dict{Symbol, Float64}=Dict{Symbol, Float64}(),
      CG::Vector{Float64}=zeros(Float64, 3))
      soln, vels=VLsolve(acft; alpha=alpha, beta=beta, p=p, q=q, r=r, control_deflections=control_deflections)

      forc=get_forces(acft, soln, vels)

      surf_results=[
        Dict(
          :sectional=>sct,
          :total=>tot
        ) for (sct, tot) in [get_surf_coefficients(surf, forc, vels; alpha=alpha, beta=beta, CG=CG) for surf in acft.surfaces]
      ]

      hinge_moments=get_hinge_moments(acft, forc)

      total_results=Dict{Symbol, Float64}()

      for (surf, sres) in zip(acft.surfaces, surf_results)
        for (k, v) in sres[:total]
          if !haskey(total_results, k)
            total_results[k]=0.0
          end

          curr=v

          if k==:Cm
            curr*=(surf.mac/acft.cref)
          elseif k==:Cn || k==:Cl
            curr*=(surf.mac/acft.bref)
          end

          curr*=(surf.S/acft.Sref)

          total_results[k]+=curr
        end
      end

      return Dict(
        :total=>total_results,
        :surface=>surf_results,
        :hinge=>hinge_moments
      )
  end
  
  export get_data

end

