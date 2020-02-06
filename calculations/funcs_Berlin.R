get_vec_prod <- function(vec_a, vec_b) 
{ 
  temp_a=dim(vec_a)
  temp_b=dim(vec_b)
  
  if(temp_a[1] != temp_b[1])
  {
   stop("vecotrs incompatible") 
  }

  if(temp_a[2] != 1 || temp_b[2]!=1)
  {stop("should be inputting column vectors")}

  result=0;    
  for(i in c(1:temp_a[1]))
  {
    result = result + vec_a[i,1]*vec_b[i,1]
  }
  
  return(result)
}



#------------------------------------
proj <- function( u, v ) #---- PROJECTS v ONTO u 
{
  temp_u = dim( as.matrix(u) )
  temp_v = dim( as.matrix(v) )
  
  if( (temp_u[2] != 1 || temp_v[2] != 1) || temp_u[1] != temp_v[1] ){print("yup")}
  
  #  { stop("input vectors incompatible with projection") }

result =  (sum(u*v) / sum(u*u)) *u
    
return(result)
}

#------------------------------------

gen_orthonorm_vecs <- function(vec_in) 
{ 
  temp=dim(vec_in)
  
  if(temp[1] <= 1 || temp[2] != 1)
  {
    stop(paste("vector input of dimension",temp[1]," x ", temp[2],"; incompatible with orthogonal generation") )
  }
  N=temp[1]
  
  result     = matrix(0,N,N)
  result[,1] = vec_in/(norm(vec_in,type="F"))
  
  for(i in c(2:N))
  {
    temp_a   = matrix(0,N,1)
    temp_b   = matrix(0,N,1)
    
    temp_a[1:(i-1)] = 1
    temp_a[i]       = 2
    
    for( j in c(1:(i-1) ) )
    {
      temp_a = temp_a - proj(result[, j, drop = FALSE],temp_a) #-- subtract projections onto each other dimension
    }
    
    result[,i] = temp_a/(norm(temp_a,type="F"))
  }
  
  return(result)
}

#============================================================================
check_Euc_dist_rows <- function(A,B) 
{
  TOL=0.001
  
  dim_A=dim(A)
  dim_B=dim(B)
  
  if(dim_A[1] != dim_B[1] || dim_A[2] != dim_B[2])
  {stop("matrices don't have conformable dimensions for Euclidean consistency check.")}
  
  num_rows=dim_A[1]
  num_ax  = dim_A[2]
  
  for(i in c(2:num_rows))
  {
    for ( j in c(1:(i-1) ) )
    {
      dist_A = norm( t(A[i,]-A[j,]), type="F")
      dist_B = norm( t(B[i,]-B[j,]), type="F")
      
      if(  (dist_A-dist_B)/dist_A > TOL )
      {      
        stop("error: distances don't agree")
      }
      else
      {
        print(paste("distance between vector ",i , " and vector ", j , "are in agreement"))
      }
    }
    
  }
}
#========================================================================================================
get_rotation_and_clusters <- function(Data,k) 
{
  dim_Data          = dim(Data)
  vote_share        = Data[,18:dim_Data[2]]/Data[,17]
  WO_array          = Data[,9];
  
  dim_vs = dim(vote_share)
  num_parties     = dim_vs[2]
  num_ridings     = dim_vs[1]

  #=======  GET THE LIST OF POINTS ON  EAST AND WEST ========
  Westlist = which(Data[,9]=="W")
  Ostlist  = which(Data[,9]=="O")
  Westpoints = vote_share[Westlist,]; names(Westpoints) = names(vote_share);
  Ostpoints  = vote_share[Ostlist,];  names(Ostpoints) = names(vote_share);
  
  #----- DO CLUSTERING ANALYSIS ----
  kmeans(vote_share,k)->kmeans_out

    if(k>2)
      { #---- the cluster number indices are random and arbitrary == shuffle them in this function so that
        #---- 1 -> "W", 2->"O", and the rest are kept as they were.
        kmeans_out = standardize_clusters(kmeans_out, k, WO_array)
      }

  frac_WO = matrix(0,k,2)
  for(i in c(1:k))
  {
    frac_WO[i,1] = sum(as.integer(as.matrix(kmeans_out$cluster)==i) * as.integer(WO_array=="W"))/(sum(as.integer(as.matrix(kmeans_out$cluster)==i) ) )
    frac_WO[i,2] = sum(as.integer(as.matrix(kmeans_out$cluster)==i) * as.integer(WO_array=="O"))/(sum(as.integer(as.matrix(kmeans_out$cluster)==i) ) )
  } #--- this should now return with the top 2x2 element subset quasi-rectangular
  
  # ---------  GET THE MEANS OF EACH -----------
  meanWest =matrix(0,1,num_parties) ; names(meanWest) = names(vote_share);
  meanOst  =matrix(0,1,num_parties) ; names(meanOst)  = names(vote_share);

  for(i in c(1:num_parties ) )
  { #--- take column by column means throughout the parties.
    meanWest [i]   = mean( Westpoints[,i]) 
    meanOst  [i]   = mean( Ostpoints[,i] )
  }

  
  if(takesep=="OW")
  {
    sep   = (meanWest-meanOst); names(sep) = names(vote_share); 
  }
  else if(takesep=="cl12")
  {
    sep =  t( kmeans_out$centers[1,]-  kmeans_out$centers[2,] )
  }
  else
  {
    stop("have not declared how to define sep")
  }
  
  #--- maybe consider redifining this as the vector seperating clusters "1" and "2"

  sepu  = sep/(sqrt(sum(sep*sep))); names(sepu) = names(sep) #--- unit vector (row) in sep direction.

  #======== ROTATE THE DATA SUCH THAT ONE BASIS VECTOR IS THE SEPERATION OF THE TWO MEANS.
  R1 = gen_orthonorm_vecs(t(sepu)) #--first rotation vector: first column runs along sepu, other columns form an arbitrary orthonormal complete basis.

  D2=t( t(R1) %*% t(vote_share) )
  # --- D2 now contains all of the points rotated in such a way that the first column is the vector pointing 
  D3      = D2[,2:num_parties] ; # Remove the first column.

  prtemp  = prcomp(D3) 
  
  D4 = t(t(prtemp$rotation) %*% t(D3))

  R2 = cbind( rbind(1,matrix(0,num_parties-1,1)) , rbind( matrix(0, 1, num_parties-1), t(prtemp$rotation)  ) ) # this is now the rotation matrix that can left multiply 
  # onto the transpose of the data matrix (i.e. column-wise points in party-space) to rotate it into the new basis.
  rot_mat = R2 %*% R1
    
  voteshare_rotated = cbind(D2[,1],D4)  #---- THIS SHOULD NOW BE OUR FINAL MATRIX. 

  # testing_dat = t( rot_mat %*% t(vote_share) )
  # check_Euc_dist_rows(testing_dat,vote_share) 
  # check_Euc_dist_rows(voteshare_rotated,vote_share)  #-- checked both of these a few times. I'm satisfied it's preserving the norm of all point-distances.
  
  Ostpoints_rotated  = voteshare_rotated [Ostlist,]
  Westpoints_rotated = voteshare_rotated [Westlist,]
  
  #===== FINISHED =======
   return( list( "rot_mat" = rot_mat, "vote_share" =  vote_share, "Datapoints_rotated" =  voteshare_rotated, "Ostpoints_rotated" =  Ostpoints_rotated, "Westpoints_rotated" = Westpoints_rotated, "sepu"=sepu, "sep"=sep, "kmeans_out"=kmeans_out) )

  }



########################################################
standardize_clusters <- function(kmeans_in, k, WO_array) 
{
  #---- the cluster number indices are random and arbitrary == shuffle them in this function so that
  #---- 1 -> "W", 2->"O", and the rest are kept as they were (if k>3).
  
  if(k<=2 || k %% 1 !=0 )
  {stop(paste("standardizing for k=", k, ",this doesnt make sense. exiting."))}

  result=kmeans_in
  
  frac_W=matrix(0,k,1)
  frac_O=matrix(0,k,1)

  #--- swap the "W" cluster with cluster "1"
  for(i in c(1:k))
  {
    frac_W[i] =sum(as.integer(as.matrix(kmeans_in$cluster)==i) * as.integer(WO_array=="W"))/(sum(as.integer(as.matrix(kmeans_in$cluster)==i) ) )
  }
  index_W = which.max(frac_W) 

  entries_1  = which(kmeans_in$cluster==1) #---the array of elements that are 1
  temp_out   = kmeans_in$cluster
  
  temp_out[which(kmeans_in$cluster==index_W)] = 1       #--- set whatever elements were labelled index_W to be labelled "1"
  temp_out[ entries_1]                        = index_W #--- set the elements that used to be 1 to index_W
  
  result$centers[1,]       = kmeans_in$centers[index_W,]
  result$centers[index_W,] = kmeans_in$centers[1,]
  
  result$withinss[1]       = kmeans_in$withinss[index_W]
  result$withinss[index_W] = kmeans_in$withinss[1]
  
  result$size[1]           =kmeans_in$size[index_W]
  result$size[index_W]     =kmeans_in$size[1]
  
  #--- Now swap the "O" cluster with cluster "2"
  for(i in c(1:k))
  {
    frac_O[i] =sum(as.integer(as.matrix(temp_out)==i) * as.integer(WO_array=="O"))/(sum(as.integer(as.matrix(temp_out)==i) ) )
  }
  index_O = which.max(frac_O) 
  
  entries_2   = which(temp_out==2) #---the array of elements that are 2 originally

  temp_out[ which(temp_out==index_O)] = 2       #--- set whatever elements were labelled index_O to be labelled "2"
  temp_out[ entries_2]                = index_O #--- set the elements that used to be 2 to the old index_O
  
  result$cluster=temp_out
  
  #---- NOW FINISH SHUFFLING THE OTHER ATTRIBUTES
  temp                      = result$centers[2,]
  result$centers[2,]        = result$centers[index_O,]
  result$centers[index_O,]  = temp
  
  temp                     = result$withinss[2]
  result$withinss[2]       = result$withinss[index_O]
  result$withinss[index_O] = temp

  temp                  = result$size[2]
  result$size[2]        = result$size[index_O]
  result$size[index_O]  = temp
  
  
  result$frac_OW=matrix(0,k,2)

  for(i in c(1:k))
    {
    result$frac_OW[i,1] = sum(as.integer(as.matrix(result$cluster)==i) * as.integer(WO_array=="W"))/(sum(as.integer(as.matrix(result$cluster)==i) ) )
    result$frac_OW[i,2] = sum(as.integer(as.matrix(result$cluster)==i) * as.integer(WO_array=="O"))/(sum(as.integer(as.matrix(result$cluster)==i) ) )
    }
  
  return(result)
}

#================================================================================================
rearrange_to_match_shape_file <- function( Data, geofile, num_Districts)
{

if( num_Districts != 12)
  { stop("unexpected number of districts")}
  
dummy         = dim(Data)
num_ridings   = dummy[1]
num_cols      = dummy[2]

used=matrix(0, num_ridings, 1)

for(i in c(1:num_Districts)) #--- we refer to the 12 larger regions as "Districts" and the smaller sublocations as "ridings"
  {
  riding_list_xl  = which(Data[,3]== i)
  riding_array_xl = as.integer( as.character( Data[riding_list_xl,5]) )
  District_Data   = Data[riding_list_xl,]
  L=length(riding_list_xl)
  
  riding_list_geo  = which(as.integer(as.character(Berlingeo_fine$BEZ)) == i)
  riding_array_geo = as.integer(as.character(Berlingeo_fine$UWB))[riding_list_geo]
  if(L!= length(riding_list_geo) || L <= 2)
    { stop("lengths of ridings from the geo shape file and vote data file dont match or are small enough to cause problems") }

    
  rot_mat_i = get_rot_mat(riding_array_geo,  riding_array_xl) # DEFINITION: first = rot_mat * second
  #--- if the resulting matrix is not orthogonal, an error flag is triggered.

  assigned        = matrix(0,L,1)
  index           = which(rot_mat_i[1,]==1)
  temp_dat        = District_Data[index,]
  assigned[index] = assigned[index] +1

  for(j in c(2:L))
    {
    index    = which(rot_mat_i[j,]==1)
    temp_dat = rbind.data.frame(temp_dat, District_Data[index,]  )
    assigned[index] = assigned[index] +1
    }
  if( min( assigned)!=1 || max(assigned)!=1)
    {stop("Error in assignment from rotation matrix; not mapping 1:1")}
  
  if(i==1)
    {
    output_mat=temp_dat  
    }
  else
    {
    output_mat = rbind.data.frame(output_mat,temp_dat)
    }
  
  used[riding_list_xl]= used[riding_list_xl]+assigned
  
  }

#result=list()
#result$output_mat = output_mat
result = output_mat
# result$used       = used

return(result)
}

#================================================================================================
district_string <- function( District )
{
  if(is.vector(District)!=1)
    stop("Need 1-D date for conversion to district strings")
  
  L=length(District)
  
  output=vector("character", L)
  
  for(i in c(1:L))
  {
    
    if(District[i] %% 1 != 0 || (District[i] <= 0 || District[i] >= 100)  )
    {stop("District[i]_string only handles integers between 1 and 100.")}
    if (District[i] >= 10 && District[i] <= 99)
    { output[i]=as.character(District[i]) }
    else if (District[i] >= 1 && District[i] <= 9)
    { output[i] = paste("0",as.character(District[i]),sep = "") }
    else
    {
      stop(paste("unexpected case: District[i] =",as.character(District[i]) ) )
    }  
    
  }
  
  
  return(output)
}  

#================================================================================================
rearrange_covariate_to_match_shape_file <- function( CD, geofile )
{
  dummy = dim(CD); num_ridings = dummy[1]; num_cols = dummy[2]
  
  compoundname_CD  = paste( district_string(CD[,2]), CD[,4], sep=""  )
  compoundname_geo = paste( geofile$BEZ, geofile$UWB, sep=""  )
  
  CD_to_geo_rot    = get_rot_mat(compoundname_geo, compoundname_CD)
  
  for(i in c(1:num_ridings)) #--- we refer to the 12 larger regions as "Districts" and the smaller sublocations as "ridings"
  {
    if(sum(CD_to_geo_rot[i,]) != 1)
      stop("rot_mat is irregular in rearranging covariate data")
    
    next_entry  = CD[ which(CD_to_geo_rot[i,]== 1), ]
    
    if(i==1)
    {
      output_mat = next_entry 
    }
    else
    {
      output_mat = rbind.data.frame(output_mat, next_entry)
    }
    
  }
  return(output_mat)
}
