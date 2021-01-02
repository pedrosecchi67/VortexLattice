module MeshPlotting

  using PyCall
  using Logging

  has_plots=false

  function __init__()
    global has_plots

    try
        py"""
        from matplotlib import pyplot as plt
        from matplotlib import cm
        from matplotlib import colors
        from mpl_toolkits import mplot3d
        import numpy as np

        def plot_surfs(pt_arrs):
            maxdim=0.0

            ax=plt.axes(projection='3d')

            for ptmat in pt_arrs:
                maxdim_loc=np.amax(np.abs(ptmat))

                if maxdim_loc>maxdim:
                    maxdim=maxdim_loc

                ax.plot_wireframe(ptmat[:, :, 0], ptmat[:, :, 1], ptmat[:, :, 2], color='b')

            ax.set_xlim((-maxdim, maxdim))
            ax.set_ylim((-maxdim, maxdim))
            ax.set_zlim((-maxdim, maxdim))

            plt.show()

        def plot_props(pt_arrs, prop_arrs):
            maxdim=0.0
            maxval=0.0
            minval=0.0

            for i, propmat in enumerate(prop_arrs):
                maxloc=np.amax(propmat)
                minloc=np.amin(propmat)

                if i==0:
                    maxval=maxloc
                    minval=minloc
                else:
                    if maxloc>maxval:
                        maxval=maxval

                    if minloc<minval:
                        minval=minloc

            ax=plt.axes(projection='3d')

            cmap=cm.get_cmap('viridis')

            for ptmat, propmat in zip(pt_arrs, prop_arrs):
                maxdim_loc=np.amax(np.abs(ptmat))

                if maxdim_loc>maxdim:
                    maxdim=maxdim_loc

                cols=np.array(
                    [
                        [list(cmap(np.interp(p, [minval, maxval], [0.0, 1.0]))) for p in line] for line in propmat
                    ]
                )

                ax.plot_surface(ptmat[:, :, 0], ptmat[:, :, 1], ptmat[:, :, 2], facecolors=cols)

            ax.set_xlim((-maxdim, maxdim))
            ax.set_ylim((-maxdim, maxdim))
            ax.set_zlim((-maxdim, maxdim))

            normalizer=colors.Normalize(vmin=minval, vmax=maxval, clip=True)
            scmap=cm.ScalarMappable(norm=normalizer, cmap=cmap)

            plt.colorbar(scmap)

            plt.show()
        """

        has_plots=true
    catch
        @warn "Unable to properly reach matplotlib via PyCall. Plotting functions will be unavailable"

        has_plots=false
    end
  end

  """
  Function to plot an aircraft's mesh geometry using matplotlib.pyplot through PyCall

  * pt_mat_arr: array of point grids
  """

  function plot_grids(pt_mat_arr::Vector{Array{Float64, 3}})
    if has_plots
        py"plot_surfs"(pt_mat_arr)
    else
        @warn "Attempting to plot when matplotlib has been detected unreachable by PyCall"
    end
end

  export plot_grids

  """
  Function to plot the values of a variable throughout the mesh

  * pt_mat_arr: array of point grids
  * props_arr: array of property values
  """

  function plot_props_in_mesh(pt_mat_arr::Vector{Array{Float64, 3}}, prop_arr::Vector{Array{Float64, 2}})
    if has_plots
        py"plot_props"(pt_mat_arr, prop_arr)
    else
        @warn "Attempting to plot when matplotlib has been detected unreachable by PyCall"
    end
  end

  export plot_props_in_mesh

end

