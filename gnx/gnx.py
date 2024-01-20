
import geonomics as gnx
import numpy as np
import pandas as pd
import multiprocessing as mp
import sys
import matplotlib.pyplot as plt
import os
from functools import partial
from os.path import exists

# Make uniform array
def make_unif_array(n):
    """Makes a square array of ones, size n x n cells."""
    array = np.ones((n, n))
    return array

def make_horizontal_split_array(n):
    """Makes a square array of size n x n cells, split horizontally.
       The top half is filled with ones, and the bottom half is filled with twos."""
    array = np.ones((n, n))  # Start with an array of ones
    half_n = n // 2  # Calculate the halfway point
    array[half_n:] = 2  # Fill the bottom half with twos
    return array

def make_vertical_split_array(n):
    """Makes a square array of size n x n cells, split vertically.
       The left half is filled with ones, and the right half is filled with 0.5."""
    array = np.ones((n, n))  # Start with an array of ones
    half_n = n // 2  # Calculate the halfway point
    array[:, half_n:] = 0.5  # Fill the right half with 0.5
    return array

klyr = make_horizontal_split_array(100)
envlyr = make_vertical_split_array(100)
unifenv = make_unif_array(100)

# Plot the arrays
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Horizontal split array
axes[0].imshow(klyr, cmap='viridis')
axes[0].set_title('Horizontal Split Array')
axes[0].axis('off')  # Turn off the axis

# Vertical split array
axes[1].imshow(envlyr, cmap='viridis')
axes[1].set_title('Vertical Split Array')
axes[1].axis('off')  # Turn off the axis

#plt.show()

params = {
    # --------------------------------------------------------------------------#

    # -----------------#
    # --- LANDSCAPE ---#
    # -----------------#
    'landscape': {

        # ------------#
        # --- main ---#
        # ------------#
        'main': {
            # x,y (a.k.a. j,i) dimensions of the Landscape
            'dim': (100, 100),
            # x,y resolution of the Landscape
            'res': (1, 1),
            # x,y coords of upper-left corner of the Landscape
            'ulc': (0, 0),
            # projection of the Landscape
            'prj': None,
        },  # <END> 'main'

        # --------------#
        # --- layers ---#
        # --------------#
        'layers': {

            # layer name (LAYER NAMES MUST BE UNIQUE!)
            'lyr_0': {

                # -------------------------------------#
                # --- layer num. 0: init parameters ---#
                # -------------------------------------#

                # initiating parameters for this layer
                'init': {

                    # parameters for a 'defined'-type Layer
                    'defined': {
                        # raster to use for the Layer
                        'rast': unifenv,
                        # point coordinates
                        'pts': None,
                        # point values
                        'vals': None,
                        # interpolation method {None, 'linear', 'cubic',
                        # 'nearest'}
                        'interp_method': None,

                    },  # <END> 'defined'

                },  # <END> 'init'

            },  # <END> layer num. 0

            # layer name (LAYER NAMES MUST BE UNIQUE!)
            'klyr': {

                # -------------------------------------#
                # --- layer num. 1: init parameters ---#
                # -------------------------------------#

                # initiating parameters for this layer
                'init': {

                    # parameters for a 'defined'-type Layer
                    'defined': {
                        # raster to use for the Layer
                        'rast': klyr,
                        # point coordinates
                        'pts': None,
                        # point values
                        'vals': None,
                        # interpolation method {None, 'linear', 'cubic',
                        # 'nearest'}
                        'interp_method': None,

                    },  # <END> 'defined'

                },  # <END> 'init'

            },  # <END> layer num. 1

            # layer name (LAYER NAMES MUST BE UNIQUE!)
            'envlyr': {

                # -------------------------------------#
                # --- layer num. 2: init parameters ---#
                # -------------------------------------#

                # initiating parameters for this layer
                'init': {

                    # parameters for a 'defined'-type Layer
                    'defined': {
                        # raster to use for the Layer
                        'rast': envlyr,
                        # point coordinates
                        'pts': None,
                        # point values
                        'vals': None,
                        # interpolation method {None, 'linear', 'cubic',
                        # 'nearest'}
                        'interp_method': None,

                    },  # <END> 'defined'

                },  # <END> 'init'

            },  # <END> layer num. 2

            #### NOTE: Individual Layers' sections can be copy-and-pasted (and
            #### assigned distinct keys and names), to create additional Layers.

        }  # <END> 'layers'

    },  # <END> 'landscape'

    # -------------------------------------------------------------------------#

    # -----------------#
    # --- COMMUNITY ---#
    # -----------------#
    'comm': {

        'species': {

            # species name (SPECIES NAMES MUST BE UNIQUE!)
            'spp_0': {

                # -----------------------------------#
                # --- spp num. 0: init parameters ---#
                # -----------------------------------#

                'init': {
                    # starting number of individs
                    'N': 10000,
                    # carrying-capacity Layer name
                    'K_layer': 'klyr',
                    # multiplicative factor for carrying-capacity layer
                    'K_factor': 2,
                },  # <END> 'init'

                # -------------------------------------#
                # --- spp num. 0: mating parameters ---#
                # -------------------------------------#

                'mating': {
                    # age(s) at sexual maturity (if tuple, female first)
                    'repro_age': 0,
                    # whether to assign sexes
                    'sex': False,
                    # ratio of males to females
                    'sex_ratio': 1 / 1,
                    # whether P(birth) should be weighted by parental dist
                    'dist_weighted_birth': False,
                    # intrinsic growth rate
                    'R': 1,
                    # intrinsic birth rate (MUST BE 0<=b<=1)
                    'b': 0.8,
                    # expectation of distr of n offspring per mating pair
                    'n_births_distr_lambda': 2,
                    # whether n births should be fixed at n_births_dist_lambda
                    'n_births_fixed': True,
                    # ADDED BY AB: choose nearest mate
                    'choose_nearest_mate': False,
                    # ADDED BY AB: choose nearest mate
                    'inverse_dist_mating': False,
                    # radius of mate-search area
                    'mating_radius': 2,
                },  # <END> 'mating'

                # ----------------------------------------#
                # --- spp num. 0: mortality parameters ---#
                # ----------------------------------------#

                'mortality': {
                    # maximum age
                    'max_age': 3,
                    # min P(death) (MUST BE 0<=d_min<=1)
                    'd_min': 0,
                    # max P(death) (MUST BE 0<=d_max<=1)
                    'd_max': 1,
                    # width of window used to estimate local pop density
                    'density_grid_window_width': None,
                },  # <END> 'mortality'

                # ---------------------------------------#
                # --- spp num. 0: movement parameters ---#
                # ---------------------------------------#

                'movement': {
                    # whether or not the species is mobile
                    'move': True,
                    # mode of distr of movement direction
                    'direction_distr_mu': 1,
                    # concentration of distr of movement direction
                    'direction_distr_kappa': 0,
                    # 1st param of distr of movement distance
                    'movement_distance_distr_param1': 0,
                    # 2nd param of distr of movement distance
                    'movement_distance_distr_param2': 0.5,
                    # movement distance distr to use
                    'movement_distance_distr': 'lognormal',
                    # 1st param of distr of dispersal distance
                    'dispersal_distance_distr_param1': 0,
                    # 2nd param of distr of dispersal distance
                    'dispersal_distance_distr_param2': 0.5,
                    # dispersal distance distr to use
                    'dispersal_distance_distr': 'lognormal',
                    'move_surf': {
                        # move-surf Layer name
                        'layer': 'lyr_0',
                        # whether to use mixture distrs
                        'mixture': True,
                        # concentration of distrs
                        'vm_distr_kappa': 12,
                        # length of approximation vectors for distrs
                        'approx_len': 5000,
                    },  # <END> 'move_surf'

                },  # <END> 'movement'

                # ---------------------------------------------------#
                # --- spp num. 0: genomic architecture parameters ---#
                # ---------------------------------------------------#

                'gen_arch': {
                    # whether to use tskit (to record full spatial pedigree)
                    'use_tskit': False,
                    # time step interval for simplication of tskit tables
                    'tskit_simp_interval': 25,  # changed from 100
                    # whether to jitter recomb bps, only needed to correctly track num_trees
                    'jitter_breakpoints': False,
                    # file defining custom genomic arch
                    # found here /p1_gnxsims/gnx/
                    'gen_arch_file': "genomic_architecture.csv",
                    # num of loci
                    'L': 10000,
                    # num of chromosomes (doesn't matter when there is no linkage)
                    'l_c': [1],
                    # starting allele frequency (None to draw freqs randomly)
                    'start_p_fixed': 0.5,
                    # whether to start neutral locus freqs at 0
                    'start_neut_zero': False,
                    # genome-wide per-base neutral mut rate (0 to disable)
                    'mu_neut': 0,
                    # genome-wide per-base deleterious mut rate (0 to disable)
                    'mu_delet': 0,
                    # shape of distr of deleterious effect sizes
                    'delet_alpha_distr_shape': 0.2,
                    # scale of distr of deleterious effect sizes
                    'delet_alpha_distr_scale': 0.2,
                    # alpha of distr of recomb rates (default = 0.5 = unlinked)
                    'r_distr_alpha': 0.5,
                    # beta of distr of recomb rates
                    'r_distr_beta': None,
                    # whether loci should be dominant (for allele '1')
                    'dom': False,
                    # whether to allow pleiotropy
                    'pleiotropy': False,
                    # custom fn for drawing recomb rates
                    'recomb_rate_custom_fn': None,
                    # number of recomb paths to hold in memory
                    'n_recomb_paths_mem': int(1e4),
                    # total number of recomb paths to simulate
                    'n_recomb_paths_tot': int(1e5),
                    # num of crossing-over events (i.e. recombs) to simulate
                    'n_recomb_sims': 10000,
                    # whether to generate recombination paths at each timestep
                    'allow_ad_hoc_recomb': False,
                    # whether to save mutation logs
                    'mut_log': False,

                    'traits': {

                        # --------------------------#
                        # --- trait 1 parameters ---#
                        # --------------------------#
                        # trait name (TRAIT NAMES MUST BE UNIQUE!)
                        'trait_1': {
                            # trait-selection Layer name
                            'layer': 'envlyr',
                            # polygenic selection coefficient
                            'phi': 1,
                            # number of loci underlying trait
                            'n_loci': 4,
                            # mutation rate at loci underlying trait
                            'mu': 0,
                            # mean of distr of effect sizes
                            'alpha_distr_mu': 0.25,
                            # variance of distr of effect size
                            'alpha_distr_sigma': 0,
                            # max allowed magnitude for an alpha value
                            'max_alpha_mag': None,
                            # curvature of fitness function
                            'gamma': 1,
                            # whether the trait is universally advantageous
                            'univ_adv': False
                        }

                        #### NOTE: Individual Traits' sections can be copy-and-pasted (and
                        #### assigned distinct keys and names), to create additional Traits.

                    },  # <END> 'traits'

                },  # <END> 'gen_arch'

            },  # <END> spp num. 0

            #### NOTE: individual Species' sections can be copy-and-pasted (and
            #### assigned distinct keys and names), to create additional Species.

        },  # <END> 'species'

    },  # <END> 'comm'

    # ------------------------------------------------------------------------#

    # -------------#
    # --- MODEL ---#
    # -------------#
    'model': {
        # total Model runtime (in timesteps)
        'T': 1001,
        # min burn-in runtime (in timesteps)
        'burn_T': 100,
        # seed number
        'num': 42,

        # -----------------------------#
        # --- iterations parameters ---#
        # -----------------------------#
        'its': {
            # num iterations
            'n_its': 1,
            # whether to randomize Landscape each iteration
            'rand_landscape': False,
            # whether to randomize Community each iteration
            'rand_comm': False,
            # whether to burn in each iteration
            'repeat_burn': False,
            #whether to randomize GenomicArchitectures each iteration
            'rand_genarch':     True,
        },  # <END> 'iterations'

        # -----------------------------------#
        # --- data-collection parameters ---#
        # -----------------------------------#
        'data': {
            'sampling': {
                # sampling scheme {'all', 'random', 'point', 'transect'}
                'scheme': 'all',
                # when to collect data
                'when': 1000,
                # whether to save current Layers when data is collected
                'include_landscape': False,
                # whether to include fixed loci in VCF files
                'include_fixed_sites': True,
            },
            'format': {
                # format for genetic data {'vcf', 'fasta'}
                'gen_format': 'vcf',
                # format for vector geodata {'csv', 'shapefile', 'geojson'}
                'geo_vect_format': 'csv',
                # format for raster geodata {'geotiff', 'txt'}
                'geo_rast_format': 'geotiff',
                #format for files containing non-neutral loci
                'nonneut_loc_format':      'csv',
            },
        },  # <END> 'data'

        # -----------------------------------#
        # --- stats-collection parameters ---#
        # -----------------------------------#
        'stats': {
            # number of individs at time t
            'Nt': {
                # whether to calculate
                'calc': True,
                # calculation frequency (in timesteps)
                'freq': 1,
            },
            # heterozgosity
            'het': {
                # whether to calculate
                'calc': True,
                # calculation frequency (in timesteps)
                'freq': 10,
                # whether to mean across sampled individs
                'mean': False,
            },
            # minor allele freq
            'maf': {
                # whether to calculate
                'calc': True,
                # calculation frequency (in timesteps)
                'freq': 10,
            },
            # mean fitness
            'mean_fit': {
                # whether to calculate
                'calc': True,
                # calculation frequency (in timesteps)
                'freq': 10,
            },
            # linkage disequilibirum
            'ld': {
                # whether to calculate
                'calc': False,
                # calculation frequency (in timesteps)
                'freq': 100,
            },
        },  # <END> 'stats'

    }  # <END> 'model'

}  # <END> params

# make our params dict into a proper Geonomics ParamsDict object
params = gnx.make_params_dict(params, "distinct_sim")
# then use it to make a model
mod = gnx.make_model(parameters=params, verbose=True)
# run the model
mod.run()