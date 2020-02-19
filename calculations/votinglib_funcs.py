import pandas as pd
from sklearn.cluster       import KMeans
from sklearn.decomposition import PCA
import math
import numpy as np

def get_major_minor_parties( DT, MAJ_PARTY_THRESH=0.05 ):
    # just using the DT (Data Table) of voting results:
    # Collect all votes
    votes_RAW       = DT.values
    n_parties_all   = DT.shape[1]
    party_names_all = DT.columns
    # Collect vote totals per party
    party_votes = [ sum( votes_RAW[:,n]) for n in range(n_parties_all) ]

    # And sum them up for the total num votes over all parties
    Total_votes_allparties = sum(party_votes)
    # print( Total_votes_allparties )

    # Now determine the fraction of all votes allocated to each party
    party_voteshare = [ party_votes[n]/Total_votes_allparties for n in range(n_parties_all) ]

    # And define the "major" parties as those above 5% total support
    major_parties = [ party_names_all[p] for p in range(n_parties_all) if party_voteshare[p]> MAJ_PARTY_THRESH ]
    minor_parties = list( party_names_all.difference( major_parties ) )

    # n_parties_major = len(major_parties)
    return (major_parties, minor_parties)

# ------------------------------------------------------
def get_normal_maj_VS( Dat_in, major_parties, minor_parties, add_other = False ):
    # Select just the major parties:
    D_major         = Dat_in[major_parties]
    D_major.index   = range( D_major.shape[0] )


    if( add_other):
        # === Collect the invalid + minor party votes, and bin them all into a separate column "other"
        Other=pd.DataFrame( { "Other": [sum((Dat_in[minor_parties+["invalid"]]).values[n,:]) for n in range(Dat_in.shape[0]) ] } )

        # Create a df with just these "parties"
        D_filtered = pd.merge( D_major,
                               Other,
                               left_index=True,
                               right_index=True )

    else:
        D_filtered = D_major

    # Normalize vote-shares
    D_normed  = D_filtered.copy()
    for i in range( D_normed.shape[0] ):
        D_normed.loc[i] = D_normed.loc[i] / sum( D_normed.loc[i] )

    return( D_normed )

# ------------------------------------------------------
def normdiff_nparrays( Dat, i, j):
    # Calculate distance between rows i, j
    return math.sqrt( sum( (Dat[i,:] - Dat[j,:])**2 ) )

# ------------------------------------------------------

def calc_percW(OWarray, kMeansarray, n):
    temp  = OWarray[ kMeansarray == n ]
    tempW = temp[temp=="W"]
    return (len(tempW)/len(temp))

def cleanclusters( k, Dat, OWarray ):
    # automate clustering, but ensure a standardized convention so that the "west" block is always 0.
    # Thereafter sorted in order of "W" percentage
    km = KMeans( n_clusters=k )
    km.fit(Dat)
    clusters = km.predict( Dat )
    # this is now an array of cluster value assignments 0..k, but the values are arbitrary and randomly assigned.
    # Now put them in order of "W" percentage so that they can be referred to consistently.

    percW     = [ calc_percW(OWarray, clusters,n) for n in set(clusters) ]
    arb_order = sorted( [ x[::-1] for x in enumerate(percW) ] )[::-1]
    # This is a set of tuples, the latter value of which indicates the arbitrary label, which we now reassign

    # Now re-order the entries in the _copy set
    clusters_ordered = clusters.copy()

    for i in range(len(percW)):
        clusters_ordered [ clusters == arb_order[i][1] ] = i

    return ( clusters_ordered, sorted( percW )[::-1] )

# ------------------------------------------------------
def iprod( u, v):
    return (sum(u*v))

# ------------------------------------------------------
def norm( u ):
    return math.sqrt( iprod(u,u) )

# ------------------------------------------------------
def dist( u, v ):
    return math.sqrt( norm(u-v))

# ------------------------------------------------------
def proj_uv( u, v ):
    #---- PROJECTS v ONTO u

    assert ( len(u.shape) == 1 and len(v.shape) ==1 )
    assert ( len(u) == len(v) )

    result =  (iprod(u,v) / norm(u)) * u

    return(result)

# ------------------------------------------------------
def gen_orthnorm_rot_mat( v ):
    assert v.shape[0] > 1 and len( v.shape )==1
    N= v.shape[0]

    result      = np.zeros([N,N])
    result[:,0] = v/(norm(v))

    # now build the rest of the vectors
    for i in range(1,N):
        temp_new    = np.ones(N)
        temp_new[i] = 2

        # ensure that they each orthogonal to all previous vectors:
        for j in range(i):
            temp_new = temp_new - proj_uv(result[:,j], temp_new)

        # And normalize them before assigning them to the output
        result[:,i] = temp_new/norm( temp_new )

    return(result)

# ========================================================
def rotate_data_along_vec_PCs( Dat, vec, TOL=0.0001 ):
    # Rotate data into, first, a preferred direction given by "vec", and then thereafter by PCs
    assert len(vec.shape) == 1
    assert Dat.shape[1]   == vec.shape [0]

    Ndim    = Dat.shape[1]
    Npoints = Dat.shape[0]

    vec_normed = (vec/norm(vec))
    Dat_output = np.zeros([Npoints, Ndim])

    # Generate a rotation matrix where the first (0) column rotates into the direction of "vec"
    # The columns after that are arbitary, so long as they all form an orthnormal complete basis
    RM0 = gen_orthnorm_rot_mat( vec_normed )

    # Now make first rotation into this basis
    Dat_rot0   = np.transpose( np.matmul( np.transpose( RM0 ), np.transpose(Dat) ) )

    # And keep only the first column
    Dat_output[:,0] = Dat_rot0[:,0]

    # The remainder will then undergo SVD:
    Dat_rot0_subset = Dat_rot0[:,1:].copy()
    U, S, V  = np.linalg.svd( Dat_rot0_subset,
                             full_matrices=False)

    # Rotate the remaining N-1 vecs along their respective PCs
    Dat_rot1 = ( V.T @ Dat_rot0_subset.T ).T
    Dat_output[:,1:] =  Dat_rot1

#    return( Dat_rot0_subset, Dat_rot1 )
    # And assign them to the remaining Dof of the output:
#    Dat_output[:,1:] = Dat_rot1
#
#     # These commented lines are sanity checks:
#     norms_0 = np.array( [ norm( Dat[n,:])   for n in range(Npoints)] )
#     norms_1 = np.array( [ norm( Dat_output[n,:]) for n in range(Npoints)] )
#     return( norms_0, norms_1)
#     plt.scatter( norms_0, norms_1 )
#     distmat_in  = np.zeros([Npoints, Npoints])
#     distmat_out = np.zeros([Npoints, Npoints])
#     for i in range(Npoints):
#        for j in range(Npoints):
#           distmat_in[i,j]  = norm( Dat[i,:]        - Dat[j,:] )
#           distmat_out[i,j] = norm( Dat_output[i,:] - Dat_output[j,:] )
#     return ( abs( distmat_out - distmat_in) )
#     # DOWN TO HERE

    return ( Dat_output )

# ========================================================
def save_fig(fig_id, tight_layout=True):
    path = os.path.join(WORKDIR, "images", fig_id + ".png")
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format='png', dpi=300)

# ========================================================
def filter_winner( Dat_table, OW_list, party_list, filterparty ):

    N=Dat_table.shape[0]

    Winner_list  = np.array( [ Dat_table.columns[np.argmax( Dat_table.iloc[n].values)] for n in  range(N) ] )

    Dat_filtered = Dat_table[ Winner_list != filterparty ]
    OW_filtered  = OW_list[   Winner_list != filterparty ]

    return( Dat_filtered, OW_filtered )

# ========================================================
def generate_ccodes( Dat_table, OW_list, party_list, pkey="GRÜNE", Wcode="blue", Ocode="red", Gcode="green"):
    # Color codes for plotting
    # takes characters or ints as code for O/W with the convention: "W"-->0, "O"-->1

    N=Dat_table.shape[0]
    # Get a list of names corresponding to the party that got the most votes
    Winner_list = np.array( [ party_list[np.argmax( Dat_table.iloc[n].values)]
                                     for n in range(N) ] )

    OW       = [ Wcode if (x=="W" or x==0 ) else  Ocode  for x in OW_list         ]
    Go20pc   = [ Gcode if (  x > 0.2      ) else "black" for x in Dat_table[pkey] ]
    Gwinner  = [ Gcode if (x==pkey or x==2) else "black" for x in Winner_list     ]

    Gwinner_OW  = [ Gcode if (Winner_list[x]==pkey or Winner_list[x]==2)
                            else Wcode if ( OW_list[x]=="W" or OW_list[x]==0 )
                            else Ocode
                    for x in range(len(Winner_list)) ]

    Go20pc_OW   = [ Gcode if Dat_table[pkey].iloc[x] >0.2
                            else Wcode if OW_list[x]=="W"
                            else Ocode
                    for x in range(len(Winner_list)) ]
    return ( {"OW"         :OW,
              "Go20pc_OW"  :Go20pc_OW,
              "Gwinner"    :Gwinner,
              "Gwinner_OW" :Gwinner_OW,
              "Go20pc_OW"  :Go20pc_OW
             } )
# ========================================================
class model_dataset():
    """A class representing a single quotient."""
    def __init__( self, Data_all, target, party_list, frac_train=0.85, rngSEED=42):

        N                = Data_all.shape[0]
        pca  = PCA()

        self.n_train     = int( frac_train * N)
        self.n_test      = N-self.n_train
        rng              = np.random.default_rng(rngSEED)
        self.train_pts   = list( rng.choice( N,
                                            size    = self.n_train,
                                            replace = False
                                           )
                               )
        self.train_dat     = Data_all.iloc[ self.train_pts ]
        self.train_dat_PCA = pca.fit_transform( self.train_dat.values )
        self.train_y       = target[ self.train_pts ]
        self.train_ccode   = generate_ccodes( self.train_dat, self.train_y, party_list)

        self.test_pts     = list ( set(range(N)).difference(set(self.train_pts)) )
        self.test_dat     = Data_all.iloc[ self.test_pts  ]
        self.test_dat_PCA = pca.fit_transform( self.test_dat.values )
        self.test_y       = target[ self.test_pts  ]
        self.test_ccode   = generate_ccodes( self.test_dat, self.test_y, party_list)

