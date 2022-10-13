


# This file contains Python3 code to accompany "Behavioral Communities and the Atomic Structure of Networks" by Matt Jackson and Evan Storms
#
# See README.txt for a description of the code
# All of the relevant python packages are included in the accompanying conda environment file 'behavioral_communities.yml'
#=============================================================================================================================================

#First import a collection of general-use python graph utility functions from the included file 'graph_utilities.py'

from graph_utilities import *


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 1: Community Detection Algorithm - Basic Function Tools 
Notes: This section defines the basic functions for testing whether a subset of a network is a convention 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#Given a partition, create a graph with the partition elements as nodes, and an edge between two elements if the elements are connected in the original graph
def partition_neighbors_graph(g,partition,partition_final):
	partition=sort_to_last(partition,partition_final)
	n=len(partition)
	m=len(partition_final)
	G=Graph(directed=False)
	G.add_vertex(n)
	for i in range(0,n-m):
		for j in range(i+1,n-m):
			for v in partition[i]:
 
				if G.edge(i,j) is None:
 
					for w in partition[j]:
 
						if g.edge(v,w) is not None:
						 
							G.add_edge(i,j)
 
							break
 
	parts=G.new_vertex_property('string')
	G.vertex_properties['contains']=parts
 
	for i in range(0,n):
		parts[i]=partition[i]
 
 
	return(G)

#Given a parition, generate a partition memebership vector for each node which keeps track of the share of the node's neighbors in each element of the partition
def gen_partition_vector(G,partition):
	partition_fractions=G.new_vertex_property('vector<double>')
	G.vertex_properties['partition_vector']=partition_fractions
	for v in G.vertices():
		partition_fractions[v]=[]
		for p in partition:
			frac = len(set.intersection(set(p),set(get_neighbors(G,v))))
			frac=frac/v.out_degree()
			partition_fractions[v].append(frac)

#Find the first element in the partition which forms a q-convention
def test_if_q_convention(G,q,partition,element=None,already_found=[]):
	#print('testing on partition', partition)
	if len(partition)==1:
		return partition[0]


	gen_partition_vector(G,partition)

 
	if element is None:
 
		for i in range(0,len(partition)) :
 
			if set(partition[i]) not in [set(a) for a in already_found]:

 
				if min([G.vp.partition_vector[j][i] for j in partition[i]]) >= q > max([G.vp.partition_vector[j][i] for j in G.vertices() if G.vertex_index[j] not in partition[i]]):
					return partition[i]

		return False
 
	else: 

		i=element

		if set(partition[i]) not in [set(a) for a in already_found]:
 
			if min([G.vp.partition_vector[j][i] for j in partition[i]]) >= q > max([G.vp.partition_vector[j][i] for j in G.vertices() if G.vertex_index[j] not in partition[i]]):
				return partition[i]
		
		return False


def check_for_q_conventions(G,q,partition,already_found=[]):
	potential_atom=test_if_q_convention(G,q,partition,already_found=already_found)
	new_already_found=deepcopy(already_found)
 
	while potential_atom is not False:

		new_already_found.append(potential_atom)

		potential_atom=test_if_q_convention(G,q,partition,already_found=new_already_found)

	return new_already_found

#Check whether the new partition contains conventions we haven't already found
def test_new_q_partition(G,q,partition,already_found,to_merge,TF=False):

	if len(to_merge) == len(partition):
		if TF is False:
			return [True,partition,already_found]
		else:
			return True 


	new_cand_partition=partition_merge_by_parts(partition,to_merge)

	new_cand_already_found=check_for_q_conventions(G,q,new_cand_partition,already_found)

	if TF is False:

		if len(new_cand_already_found)>len(already_found):
			#print('updating partition')
			partition=new_cand_partition
			already_found=new_cand_already_found
	
			partition=sort_to_last(partition,already_found)


			return [True,partition,already_found]

		return [False,partition,already_found]

	else:

		if len(new_cand_already_found)>len(already_found):
			#print('updating partition')
			partition=new_cand_partition
			already_found=new_cand_already_found
	
			partition=sort_to_last(partition,already_found)


			return True

		return False




""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 2: Functions for finding the atomic structure of a graph by growing small subsets 

Notes: This section implements the approximation algorithm described in the text for finding the q-atoms of a network by growing 
small subsets of the network into conventions. 

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""








#Find minimal q-closed superset of S 
def find_minimal_q_closed_superset(S,G,q):
	S=list(S)
	S_complement=node_complements(G,S)
	added=False 

	for v_num in S_complement:
		v=G.vertex(v_num)
		v_share_in_S=len(set.intersection(set(S),set(get_neighbors(G,v))))/v.out_degree()
		if v_share_in_S>=q:
			S = S + [v_num]
			added=True

	if added is True:
		return find_minimal_q_closed_superset(S,G,q)

	else:
		return frozenset(S) 



#Take the q closure of each set in collection
def take_q_closures(collection,G,q,num_cores=30,clear_main=False,seed_dict=None,PBT=None,PB=None,MPB=None,proc_name='Generating q closures'):
	find_minimal_q_closed_superset_mod=partial(find_minimal_q_closed_superset,G=G,q=q)
	if seed_dict is None:
		results=many_pool_imap_mini(find_minimal_q_closed_superset_mod,collection,proc_name=proc_name,num_cores=num_cores,clear_main=clear_main,PB=PB,MPB=MPB,PBT=PBT)
		return results 
	else:
		results,new_dict=many_pool_imap_mini_w_dict(find_minimal_q_closed_superset_mod,collection,origin_dict=seed_dict,proc_name=proc_name,num_cores=num_cores,clear_main=clear_main,PB=PB,MPB=MPB,PBT=PBT)
		return results,new_dict 


#Add nodes to a subset S until it is q-cohesive (works by adding the nodes which most increases the subset's cohesion, taking the q-closure, then iterating until we get a cohesive set)
def cohesify(S,G,q):

	if len(S)==0:
		return None

	if len(S)==G.num_vertices():
		return S

	g=G.copy()
	S_complement=node_complements(G,S)
	gen_partition_vector(g,[S,S_complement])
	min_cohesion=min([g.vp.partition_vector[v][0] for v in S])
	#print('min cohesion is',min_cohesion)

	if min_cohesion<q:
		u= min([v for v in S if g.vp.partition_vector[v][0] == min_cohesion])
		u_neighbours= [g.vertex_index[v] for v in list(g.vertex(u).out_neighbours()) if g.vertex_index[v] not in S]
		max_into = max([g.vp.partition_vector[v][0] for v in u_neighbours])
		potential_to_adds=[v for v in u_neighbours if g.vp.partition_vector[v][0] ==max_into]
		#print('potential new additions are',potential_to_adds)

		#If there are multiple potential additions, take the one which gives the smallest closure 
		if len(potential_to_adds)>0:
			closures=[find_minimal_q_closed_superset(set(list(S) + [a]),G,q) for a in potential_to_adds]
			nodes_and_closures=list(zip(potential_to_adds,closures))
			min_closure_size=min([len(x) for x in closures])
			keep=[a for a in nodes_and_closures if len(a[1])==min_closure_size]
			nodes=[a[0] for a in keep]
			if len(nodes)>1:
				cohesions=[find_subset_cohesion(a[1],g) for a in keep]
				nodes_and_cohesions=zip(nodes,cohesions)
				max_cohesion=max(cohesions)
				keep=[a for a in nodes_and_cohesions if a[1]==max_cohesion]
				nodes=[a[0] for a in keep]
			to_add=min(nodes)
		else:
			to_add =potential_to_adds[0]

		S = set(list(S) + [to_add])
		#print('new S is',S)
		return cohesify(find_minimal_q_closed_superset(S,G,q),G,q)

	else: 
		return S 


#Function that returns the cohesion (minimal in-subset share of neighbors) for a subset S
def find_subset_cohesion(S,G):
	g=G.copy()
	S_complement=node_complements(G,S)
	gen_partition_vector(g,[S,S_complement])
	min_cohesion=min([g.vp.partition_vector[v][0] for v in S])
	return min_cohesion



#Test whether a subset is q-cohesive
def test_cohesion(S,G,q):
	if len(S)==0:
		return None 

	g=G.copy()
	S_complement=list(node_complements(g,S))
	gen_partition_vector(g,[S,S_complement])
	min_cohesion=min([g.vp.partition_vector[v][0] for v in S])

	if min_cohesion>=q:
		return True
	else:
		return False 


#Test whether a subset is q-closed 
def test_closed(S,G,q):

	if len(S)==0:
		return None

	if len(S) == G.num_vertices():
		return True
	else: 
		g=G.copy()
		S_complement=list(node_complements(g,S))
		gen_partition_vector(g,[S,S_complement])
		max_in=max([g.vp.partition_vector[v][0] for v in S_complement])
		if max_in < q:
			return True
		else:
			return False



#Test whether a subset is a q-convention 
def test_convention(S,G,q):
	if len(S)==0:
		return None 

	if test_cohesion(S,G,q) and test_closed(S,G,q):
		return True
	else:
		return False




#Given a subset, grow it into a q-convention
def conventionify(S,G,q,already_closed=False):
	if already_closed is False:
		S=find_minimal_q_closed_superset(S,G,q)
	cohesive=False
	while cohesive is False:
		S=find_minimal_q_closed_superset(cohesify(S,G,q),G,q)
		cohesive = test_cohesion(S,G,q)
	return S



#Return the minimal q-closures of each connected set of nodes of size less than or equal to k
def gen_minimal_q_closures_of_connected_subsets(G,q,k,blocks=[]):
	lower_admission=floor(q * min([v.out_degree() for v in G.vertices()]))

	blocked_nodes=set([])
	for block in blocks:
		blocked_nodes=blocked_nodes.union(block)


	partition=blocks + [set([G.vertex_index[v]]) for v in G.vertices() if G.vertex_index[v] not in blocked_nodes ]
	pgraph=partition_neighbors_graph(G,partition,[])
	minimal_q_closeds=[]
	subsets_pre_translate = find_k_connected_subsets_v_by_v(pgraph,k)
	subsets=[translate([a],partition)[0] for a in subsets_pre_translate]
	minimal_q_closeds=[a for a in subsets if len(a) < lower_admission]
	to_check=[a for a in subsets if len(a) >= lower_admission]
	for subset in to_check:
		candidate=find_minimal_q_closed_superset(subset,G,q)
		if candidate not in minimal_q_closeds:
			minimal_q_closeds.append(candidate)

	return list(convert_to_set_of_sets(minimal_q_closeds))




#To save time, we start by generating pre-atoms, where two nodes are the in the same pre-atom if they cannot be separated by the minimal q-closed superset of some two-node set 

def gen_q_pre_atoms(G,q):
	minimal_q_closeds=gen_minimal_q_closures_of_connected_subsets(G,q,2)
	complements=[]
	for subset in minimal_q_closeds:
		complements.append(node_complements(G,subset))

	cuts=minimal_q_closeds + complements

	pre_atoms=get_intersections(G,cuts)

	return pre_atoms



#Generate q-conventions by growing all connected subsets of size less than or equal to k 
def gen_some_q_conventions(G,q,k,blocks=[]):
	q_closeds=gen_minimal_q_closures_of_connected_subsets(G,q,k,blocks)
	conventions=[]
	for q_closed in q_closeds: 
		candidate= conventionify(q_closed,G,q,already_closed=True)
		if candidate not in conventions:
			conventions.append(candidate)

	return conventions 


#Like the above, but this stops looking for conventions once we have enough conventions to generate the singleton atomic structure 
def gen_enough_q_conventions(G,q,k,blocks=[]):
	q_closeds=gen_minimal_q_closures_of_connected_subsets(G,q,k,blocks)
	conventions=[]
	n=G.num_vertices()
	for q_closed in q_closeds: 
		candidate= conventionify(q_closed,G,q,already_closed=True)
		if candidate not in conventions:
			conventions.append(candidate)
			if len(conventions)%(2*n)==0:
				complements=[node_complements(G,a) for a in conventions]
				cuts = conventions + complements 
				atoms=get_intersections(G,cuts)
				if len(atoms)==n:
					print('found shortcut')
					return conventions

	return conventions 



#See above
def gen_some_q_conventions_multicore(G,q,k,num_cores=30,blocks=[]):
	q_closeds=gen_minimal_q_closures_of_connected_subsets_multicore(G,q,k,blocks)
	conventions=[]
	n=G.num_vertices()
	pool=mp.Pool(num_cores)
	conventionify_mod=partial(conventionify,G=G,q=q,already_closed=True)
	new_conventions=many_pool_imap(conventionify_mod,q_closeds,proc_name='Converting q-closed sets to conventions')
	conventions=list(convert_to_set_of_sets(new_conventions))

	return conventions 


#Generate the approximate q-atomic structure by algebra generated by the conventions obtained by growing all the subsets of size k or less into conventions (shortcut=this version stops once there are enough conventions to get complete fragmentaiton)
def gen_q_atoms_approx_shortcut(G,q,k):
	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	nulls=[G.vertex_index[v] for v in G.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return nulls
	else:
		pre_atoms=gen_q_pre_atoms(H,q)
		conventions=gen_enough_q_conventions(H,q,k,pre_atoms)
		complements=[node_complements(H,a) for a in conventions]
		cuts = conventions + complements + [{a} for a in nulls]
		return get_intersections(G,cuts)




def gen_q_conventions_approx(G,q,k):
	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	nulls=[G.vertex_index[v] for v in G.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return nulls
	else:
		pre_atoms=gen_q_pre_atoms(H,q)
		conventions=gen_some_q_conventions(H,q,k,pre_atoms)
		return conventions




def gen_q_atoms_approx(G,q,k):
	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	nulls=[G.vertex_index[v] for v in G.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return nulls
	else:
		pre_atoms=gen_q_pre_atoms(H,q)
		conventions=gen_some_q_conventions(H,q,k,pre_atoms)
		complements=[node_complements(H,a) for a in conventions]
		cuts = conventions + complements +  [{a} for a in nulls]
		return get_intersections(G,cuts)



"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 3: Generate the atomic structure by brute force checking all connected subsets for conventions. 
Notes: Not feasible for n>30, this is really just a check that the above algorithms work. 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


def gen_conventions(G,blocks=[],max_merge=None):
	blocked_nodes=[]
	for block in blocks:
		blocked_nodes=blocked_nodes + block 

	partition=blocks + [[G.vertex_index[v]] for v in G.vertices() if G.vertex_index[v] not in blocked_nodes ]
	conventions=[[i] for i in range(0,len(blocks))]

	if max_merge is None:
		max_merge=min(len(partition),15)
	current_merge=1
	last_collection=[set([i]) for i in range(0,len(partition))]
	pgraph=partition_neighbors_graph(G,partition,[])

	while current_merge<max_merge:
		last_collection=find_plus_one_connected_subsets(pgraph,last_collection)

		for S in last_collection:
			#print('testing', S)
			if test_new_partition(G,partition,[],blocks,S,TF=True):
				conventions.append(S) 

		current_merge=current_merge + 1 
		print('moving to merges of size ', current_merge)

	conventions=translate(conventions,partition)

	if set([G.vertex_index[v] for v in G.vertices()]) not in conventions:
		conventions.append(set([G.vertex_index[v] for v in G.vertices()]))

	return conventions





def gen_q_conventions(G,q,blocks=[],max_merge=None):

	blocked_nodes=[]
	for block in blocks:
		blocked_nodes=blocked_nodes + block 


	partition=blocks + [[G.vertex_index[v]] for v in G.vertices() if G.vertex_index[v] not in blocked_nodes ]
	conventions=[[i] for i in range(0,len(blocks))]

	if max_merge is None:
		max_merge=min(len(partition),15)

	current_merge=1
	last_collection=[set([i]) for i in range(0,len(partition))]
	pgraph=partition_neighbors_graph(G,partition,[])

	while current_merge<max_merge:
		last_collection=find_plus_one_connected_subsets(pgraph,last_collection)

		for S in last_collection:
			if len(S)>=len(partition):
				conventions.append(S)
				break 
			#print('testing', S)
			if test_new_q_partition(G,q,partition,blocks,S,TF=True):
				conventions.append(S) 
				#print('adding',S)
			

		current_merge=current_merge + 1 
		print('moving to merges of size ', current_merge)

	conventions=translate(conventions,partition)

	if set([G.vertex_index[v] for v in G.vertices()]) not in conventions:
		conventions.append(set([G.vertex_index[v] for v in G.vertices()]))

	return conventions


def gen_atoms(G,blocks=[],max_merge=None):

	conventions=gen_conventions(G,blocks,max_merge)


	return get_intersections(G,conventions)


def gen_q_atoms(G,q,blocks=[],max_merge=None):

	conventions=gen_q_conventions(G,q,blocks,max_merge)
	complements=[]
	for convention in conventions:
		complements.append([G.vertex_index[v] for v in G.vertices() if G.vertex_index[v] not in convention])

	conventions=conventions + complements

	return get_intersections(G,conventions)



"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 4: Community Detection Algorithm - q Ranges
Notes: This section adopts the algorithm of Section 3 to find q-range robust conventions. Approach is to basically do the 
same as before at the low q, but only keep a convention if it is also cohesive at the high q. 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#here q range is a list [lower_bound,upper_bound]

def gen_some_q_range_conventions(G,q_range,k,blocks=[]):
	q_lower=q_range[0]
	q_upper=q_range[1]
	q_closeds=gen_minimal_q_closures_of_connected_subsets(G,q_lower,k,blocks)
	conventions=[]
	for q_closed in q_closeds: 
		candidate= conventionify(q_closed,G,q_lower,already_closed=True)
		if candidate not in conventions and candidate is not None and test_cohesion(candidate,G,q_upper) is True:
			conventions.append(candidate)

	return conventions 


def gen_enough_q_range_conventions(G,q_range,k,blocks=[]):
	n=G.num_vertices()
	q_lower=q_range[0]
	q_upper=q_range[1]
	q_closeds=gen_minimal_q_closures_of_connected_subsets(G,q_lower,k,blocks)
	conventions=[]
	for q_closed in q_closeds: 
		candidate= conventionify(q_closed,G,q_lower,already_closed=True)
		if candidate not in conventions and candidate is not None and test_cohesion(candidate,G,q_upper) is True:
			conventions.append(candidate)
			if len(conventions)%(2*n)==0:
				complements=[node_complements(G,a) for a in conventions]
				cuts = conventions + complements 
				atoms=get_intersections(G,cuts)
				if len(atoms)==n:
					print('found shortcut')
					return conventions

	return conventions 



def gen_q_range_atoms_approx(G,q_range,k):

	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	nulls=[G.vertex_index[v] for v in G.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return nulls
	else:
		q_lower=q_range[0]
		q_upper=q_range[1] 
		pre_atoms=gen_q_pre_atoms(H,q_lower)
		conventions=gen_some_q_range_conventions(H,q_range,k,pre_atoms)
		complements=[node_complements(H,a) for a in conventions]
		cuts = conventions + complements + [[a] for a in nulls]
		return get_intersections(G,cuts)





def gen_q_range_atoms_approx_multicore(G,q_range,k):
	PBT=task_PB(proc_name='Finding the atoms for q_range= ' + format_number(q_range[0]) + '   to  ' + format_number(q_range[1]))
	PBT.setup()
	clear_mega_progress()
	clear_main_progress()
	clear_mini_progress()

	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	nulls=[G.vertex_index[v] for v in G.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return nulls
	else:
		q_lower=q_range[0]
		q_upper=q_range[1] 
		pre_atoms=gen_q_pre_atoms_multicore(H,q_lower)
		PBT.update_runtime()
		conventions=gen_some_q_range_conventions_multicore(H,q_range,k,pre_atoms)
		PBT.update_runtime()

		complements=[node_complements(H,a) for a in conventions]
		cuts = list(conventions) + complements + [[a] for a in nulls]
		PBT.update_runtime()

		return get_intersections(G,cuts)





def gen_q_range_atoms_approx_shortcut(G,q_range,k):
	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	nulls=[G.vertex_index[v] for v in G.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return nulls
	else:
		q_lower=q_range[0]
		q_upper=q_range[1] 
		pre_atoms=gen_q_pre_atoms(H,q_lower)
		conventions=gen_enough_q_range_conventions(H,q_range,k,pre_atoms)
		complements=[node_complements(H,a) for a in conventions]
		cuts = conventions + complements + [[a] for a in nulls]
		return get_intersections(G,cuts)




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 5: Random Convention Generating Procedures 
Notes: The section gives some functions for randomly generating a q-threshold coordination eq (possibly noised up) on G. 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def seed_randomly(G,alpha):
	on=[]

	if alpha<1:
		for v in G.vertices():
			if random.random()<alpha:
				on.append(G.vertex_index[v])
	else:
		on= random.sample([G.vertex_index[v] for v in G.vertices()],alpha)
	return on



def update_behavior(G,q,on):
	new_on=[]
	for v in G.vertices():
		N_v=[G.vertex_index[u] for u in G.vertex(v).out_neighbours()]
		share_in=len([a for a in N_v if a in on])/len(N_v)
		if share_in >= q:
			new_on.append(G.vertex_index[v])

	return new_on 



def iteratively_update(on,G,q,permanent=True,max_iters=1000):
	iters=0
	on=set(on)
	if permanent is False: 
		while test_closed(on,G,q) is False and iters<max_iters:
			new_on=update_behavior(G,q,on)
			if set(on)==set(new_on):
				break
			else:
				on=new_on
			iters=iters + 1 
	else:
		return find_minimal_q_closed_superset(on,G,q) 

	if iters==max_iters:
		return None
	else:
		return on

def iteratively_update_return_set(on,G,q,permanent=True,max_iters=1000):
	out=iteratively_update(on,G,q,permanent=permanent,max_iters=max_iters)
	if out is not None:
		return frozenset(out)
	else:
		return(out)


def iteratively_update_to_find_convention(on,G,q,permanent=False,max_iters=1000):
	iters=0
	if permanent is False: 
		while test_convention(on,G,q) is False and iters<max_iters:
			new_on=update_behavior(G,q,on)
			iters=iters + 1
			if set(on)==set(new_on):
				break
			else:
				on=new_on 
	else:
		while test_convention(on,G,q) is False and iters<max_iters:
			new_on=list(set(on=update_behavior(G,q,on) + on ))
			iters=iters + 1
			if set(on)==set(new_on):
				break
			else:
				on=new_on 

	if iters==max_iters:
		return None
	else:
		return on




def generate_random_conventions(G,q,alpha):
	on=seed_randomly(G,alpha)
	return iteratively_update_to_find_convention(on,G,q)


def generate_random_conventions_just_grow(G,q,k):
	on=seed_randomly(G,alpha=k)
	return conventionify(on,G,q)





 
def generate_conventions_with_noise(G,q,alpha,beta,start_prob):
	independents,independents_on=get_independents(G,alpha,beta)
	on=seed_randomly_with_noise(G,start_prob,independents)

	last_on=[9999999999999999999]
	iters=100
	counter=0

	while last_on != on and counter<iters:  
		last_on=on
		on=update_behavior_with_noise(G,q,on,independents,independents_on)
		counter=counter + 1 

	if counter==iters:
		return []
	else:
		return on+independents_on 




def noise_up_from_convention(G,q,alpha,beta,num_to_find=4):
	starting_convention=gen_a_particular_sized_q_convention(G,q,2,1/10,2/5)
	if len(starting_convention)==0:
		return []
	independents,independents_on=get_independents(G,alpha,beta)

	on=[a for a in starting_convention if a not in independents] 
	last_on=[9999999999999999999]
	iters=100
	counter=0

	while last_on != on and counter<iters:  
		last_on=on
		on=update_behavior_with_noise(G,q,on,independents,independents_on)
		counter=counter + 1 

	if counter==iters:
		return []
	else:
		return on+independents_on 




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 6: Optimal, random, and atom-based seeding procedures
Notes: The section gives the functions which implement optimal, random, and atob-based seeding procedures (i.e. the atom-based heuristic)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



#Test one instance of random seeding in a network G (here alpha=number seeds)
def random_seed_result(G,alpha,q,permanent=True):
	on=seed_randomly(G,alpha)
	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	isoalted_seeds=[a for a in on if has_neighbor[a]==0]
	on_H=[a for a in on if has_neighbor[a]==1]

	reach=iteratively_update(on_H,H,q,permanent=permanent)
	return len(reach)+len(isoalted_seeds)



#Find the optimal seeding by brute force
def optimal_seeding_by_brute_force(G,q,k,permanent=True,start_from=None):
	vertices_to_consider=[G.vertex_index[v] for v in G.vertices() if G.vertex(v).out_degree()>ceil(1/q)]

	if start_from is None:
		combs=[list(a) for a in list(combinations(vertices_to_consider,k))]
	else:
		combs=start_from
		print("on main call avoided doing ", nCr(len(vertices_to_consider),k) - len(combs), " redundant checks") 

	starting_length=len(combs)

	max_found=0 
	current_maximizer=[]
	while len(combs)>0 and max_found<G.num_vertices():
		comb=combs.pop(0)
		reach=len(iteratively_update(comb,G,q,permanent=permanent))
		if reach>max_found:
			max_found=reach
			current_maximizer=comb

	return [max_found,current_maximizer]




#The remaining functions in this section are used in the atom-based seeding heurisitc


def pick_hard_stop(l,k):
	while nCr(l,k)>10**6:
		k=k-1
	return k 



def seeding_cost(S,G,q,num_seeds,permanent=True,start_at=1):
	j=start_at
	upper_bound=min(len(S),num_seeds+1)
	if len(S)==1:
		return [1,list(S)]
	while j<upper_bound:
		for comb in combinations(S,j):
			final_ons=iteratively_update(list(comb),G,q,permanent=permanent)
			if set(final_ons).intersection(set(S)) ==set(S):
				return [j,list(comb)]
		j=j+1 


	if j==len(S):
		return [j,list(S)] 

	if j==num_seeds+1:
		return [num_seeds +1, list(S)] 




def seed_atoms_greedily(subsets,cost_dict,k):
	seed_set=[]
	working_subsets=deepcopy(subsets)
	while len(seed_set)<k and len(working_subsets)>0:
		consider=working_subsets.pop(0)
		if cost_dict[consider][0]<= k -  len(seed_set):
			seed_set = seed_set + list(cost_dict[consider][1])

	return seed_set
 



#If we don't use all k seeds, just select k nodes that aren't in the reach
def find_leftover_seeds(G,q,seed_set,k,permanent=True):
	has_neighbor=nodes_with_neighbor(G)
	H=GraphView(G,vfilt=has_neighbor)

	already_hit=iteratively_update(seed_set,H,q,permanent=permanent)
	if len(already_hit)==G.num_vertices():
		return []

	num_leftover_seeds=k - len(seed_set)
	leftover_seeds=list(node_complements(G,already_hit))[0:num_leftover_seeds]
	return leftover_seeds 


def atom_based_heuristic(G,q,k,permanent=True):
	atoms=gen_q_atoms_approx_shortcut(G,q,3)
	costs={}

	for atom in atoms:
		if len(atom)>1:
			costs[atom]=seeding_cost(atom,G,q,permanent=permanent)
		else:
			costs[atom]=1

	atoms=sorted(atoms,key= lambda x: -len(x)/costs[x][0])
	seed_set=seed_atoms_greedily(atoms,costs,k)
	if len(seed_set)<k:
		seed_set = seed_set + find_leftover_seeds(G,q,seed_set,k,permanent=permanent)

	return [len(iteratively_update(seed_set,G,q,permanent=permanent)),seed_set]
 





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 7:T Version of Algorithm
Notes: Everything here is teh parallel of the q case in Section 3, just for the absolute threshold T. Only real difference
is that we first have to the prune the network to the T-core, since everything outside the T-core is always off. 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 7-1 T Version: Community Detection Algorithm - Basic Function Tools (Absolute Threshold Version)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

def prune_T(G,T):
	g=G.copy()
	min_degree=min([v.out_degree() for v in g.vertices()])
	while min_degree<T:
		nodes_to_be_removed=[g.vertex_index[v] for v in g.vertices() if v.out_degree() < T]
		print('nodes to remove are', nodes_to_be_removed)
		g=remove_all_edges_from_nodes(g,nodes_to_be_removed)
		if g.num_edges()==0:
			return g 
		else:
			min_degree=min([v.out_degree() for v in g.vertices() if v.out_degree()>0])
			print(min_degree)

	return g 



def gen_partition_vector_T(G,partition):
	partition_numbers=G.new_vertex_property('vector<double>')
	G.vertex_properties['partition_vector']=partition_numbers
	for v in G.vertices():
		partition_numbers[v]=[]
		for p in partition:
			num = len(set.intersection(set(list(p)),set(get_neighbors(G,v))))
			partition_numbers[v].append(num)



def test_if_T_convention(G,T,partition,element=None,already_found=[]):
	#print('testing on partition', partition)
	if len(partition)==1:
		return partition[0]


	gen_partition_vector_T(G,partition)

 
	if element is None:
 
		for i in range(0,len(partition)) :
 
			if set(partition[i]) not in [set(a) for a in already_found]:

 
				if min([G.vp.partition_vector[j][i] for j in partition[i]]) >= T > max([G.vp.partition_vector[j][i] for j in G.vertices() if G.vertex_index[j] not in partition[i]]):
					return partition[i]

		return False
 
	else: 

		i=element

		if set(partition[i]) not in [set(a) for a in already_found]:
 
			if min([G.vp.partition_vector[j][i] for j in partition[i]]) >= T > max([G.vp.partition_vector[j][i] for j in G.vertices() if G.vertex_index[j] not in partition[i]]):
				return partition[i]
		
		return False



def check_for_T_conventions(G,T,partition,already_found=[]):
	potential_atom=test_if_T_convention(G,T,partition,already_found=already_found)
	new_already_found=deepcopy(already_found)
 
	while potential_atom is not False:

		new_already_found.append(potential_atom)

		potential_atom=test_if_T_convention(G,T,partition,already_found=new_already_found)

	return new_already_found



def test_new_T_partition(G,T,partition,already_found,to_merge,TF=False):

	if len(to_merge) == len(partition):
		if TF is False:
			return [True,partition,already_found]
		else:
			return True 


	new_cand_partition=partition_merge_by_parts(partition,to_merge)

	new_cand_already_found=check_for_T_conventions(G,T,new_cand_partition,already_found)

	if TF is False:

		if len(new_cand_already_found)>len(already_found):
			#print('updating partition')
			partition=new_cand_partition
			already_found=new_cand_already_found
	
			partition=sort_to_last(partition,already_found)


			return [True,partition,already_found]

		return [False,partition,already_found]

	else:

		if len(new_cand_already_found)>len(already_found):
			#print('updating partition')
			partition=new_cand_partition
			already_found=new_cand_already_found
	
			partition=sort_to_last(partition,already_found)


			return True

		return False




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 7-2 T Version: Community Detection Algorithm
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



#Find minimal T-closed superset of S 
def find_minimal_T_closed_superset(S,G,T):
	S=list(S)
	g=G.copy()
	S_complement=node_complements(g,S)
	gen_partition_vector_T(g,[S,S_complement])

	added=False

	for v in S_complement:
		if g.vp.partition_vector[v][0]>=T:
			S = S + [g.vertex_index[v]]
			added=True

	if added is True:
		#print('New call, and S is', S)
		return find_minimal_T_closed_superset(S,G,T)

	return set(S)


def cohesify_T(S,G,T):

	if len(S)==0:
		return None

	if len(S)==G.num_vertices():
		return S

	g=G.copy()
	S_complement=node_complements(G,S)
	gen_partition_vector_T(g,[S,S_complement])
	min_cohesion=min([g.vp.partition_vector[v][0] for v in S])
	#print('min cohesion is',min_cohesion)

	if min_cohesion<T:
		u= min([v for v in S if g.vp.partition_vector[v][0] == min_cohesion])
		u_neighbours= [g.vertex_index[v] for v in list(g.vertex(u).out_neighbours()) if g.vertex_index[v] not in S]
		max_into = max([g.vp.partition_vector[v][0] for v in u_neighbours])
		potential_to_adds=[v for v in u_neighbours if g.vp.partition_vector[v][0] ==max_into]
		#print('potential new additions are',potential_to_adds)

		#If there are multiple potential additions, cut down by smallest closure then most cohesive closure
		if len(potential_to_adds)>1:
			closures=[find_minimal_T_closed_superset(set(list(S) + [a]),G,T) for a in potential_to_adds]
			nodes_and_closures=list(zip(potential_to_adds,closures))
			min_closure_size=min([len(x) for x in closures])
			keep=[a for a in nodes_and_closures if len(a[1])==min_closure_size]
			nodes=[a[0] for a in keep]
			if len(nodes)>1:
				cohesions=[find_subset_cohesion_T(a[1],g) for a in keep]
				nodes_and_cohesions=zip(nodes,cohesions)
				max_cohesion=max(cohesions)
				keep=[a for a in nodes_and_cohesions if a[1]==max_cohesion]
				nodes=[a[0] for a in keep]
			to_add=min(nodes)
		else:
			to_add =potential_to_adds[0]

		S = set(list(S) + [to_add])
		#print('new S is',S)
		return cohesify_T(find_minimal_T_closed_superset(S,G,T),G,T)

	else: 
		return S 






def find_subset_cohesion_T(S,G):
	g=G.copy()
	S_complement=[g.vertex_index[v] for v in g.vertices() if g.vertex_index[v] not in S]
	gen_partition_vector_T(g,[S,S_complement])
	min_cohesion=min([g.vp.partition_vector[v][0] for v in S])

	return min_cohesion


def test_cohesion_T(S,G,T):
	min_cohesion=find_subset_cohesion_T(S,G)
	if min_cohesion>=T:
		return True
	else:
		return False 


def test_closed_T(S,G,T):
	if len(S)==0:
		return None

	if len(S) == G.num_vertices():
		return True
	else: 
		g=G.copy()
		S_complement=list(node_complements(g,S))
		gen_partition_vector_T(g,[S,S_complement])
		max_in=max([g.vp.partition_vector[v][0] for v in S_complement])
		if max_in < T:
			return True
		else:
			return False

def test_convention_T(S,G,T):
	if len(S)==0:
		return None 

	if test_cohesion_T(S,G,T) and test_closed_T(S,G,T):
		return True
	else:
		return False









def conventionify_T(S,G,T,already_closed=False):	
	#print('starint S is',S)
	if already_closed is False:
		S=find_minimal_T_closed_superset(S,G,T)
		#print('now S is',S)

	cohesive=False
	while cohesive is False:
		S=find_minimal_T_closed_superset(cohesify_T(S,G,T),G,T)
		#print('now S is',S)

		cohesive = test_cohesion_T(S,G,T)

	return S




def gen_minimal_T_closeds_of_pairs(G,T):
	minimal_T_closeds=[]
	for edge in G.edges():
		S=[G.vertex_index[edge.source()],G.vertex_index[edge.target()]]

		candidate=find_minimal_T_closed_superset(S,G,T)
		if candidate not in minimal_T_closeds:
			minimal_T_closeds.append(candidate)

	return minimal_T_closeds



def gen_minimal_T_closeds_of_subsets(G,T,k,blocks=[]):

	if len(blocks)>0:
		blocked_nodes=set([])
		for block in blocks:
			blocked_nodes=blocked_nodes.union(block)

		partition=blocks + [set([G.vertex_index[v]]) for v in G.vertices() if G.vertex_index[v] not in blocked_nodes ]


	need_to_close=[a for a in subsets if len(a)>=T]

	minimal_T_closeds=already_in
	for S in need_to_close:
		closure=find_minimal_T_closed_superset(S,G,T)
		if closure not in minimal_T_closeds:
			minimal_T_closeds.append(closure)


	return minimal_T_closeds


def gen_some_T_conventions(G,T,k,blocks=[]):
	T_closeds=gen_minimal_T_closeds_of_subsets(G,T,k,blocks)
	conventions=[]
	for T_closed in T_closeds: 
		candidate= conventionify_T(T_closed,G,T,already_closed=True)
		if candidate not in conventions:
			conventions.append(candidate)

	return conventions 


def gen_enough_T_conventions(G,T,k,blocks=[]):
	T_closeds=gen_minimal_T_closeds_of_subsets(G,T,k,blocks)
	conventions=[]
	n=G.num_vertices()
	for T_closed in T_closeds: 
		candidate= conventionify_T(T_closed,G,T,already_closed=True)
		if candidate not in conventions:
			conventions.append(candidate)
			if len(conventions)%(2*n)==0:
				complements=[node_complements(G,a) for a in conventions]
				cuts = conventions + complements 
				atoms=get_intersections(G,cuts)
				if len(atoms)==n:
					print('found shortcut')
					return conventions

	return conventions 




def gen_T_atoms_approx_shortcut(G,T,k):
	G_prime=prune_T(G,T)
	has_neighbor=nodes_with_neighbor(G_prime)
	H=GraphView(G_prime,vfilt=has_neighbor)

	nulls=[G_prime.vertex_index[v] for v in G_prime.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return [nulls]
	else:	
		conventions=gen_enough_T_conventions(H,T,k,[])
		complements=[node_complements(H,a) for a in conventions]
		cuts = conventions + complements + [nulls]
		return get_intersections(G,cuts)

def gen_T_conventions_approx(G,T,k):
	H=prune_T(G,T)
	conventions=gen_some_T_conventions(H,T,k,[])
	return conventions




def gen_T_atoms_approx(G,T,k):

	G_prime=prune_T(G,T)
	has_neighbor=nodes_with_neighbor(G_prime)
	H=GraphView(G_prime,vfilt=has_neighbor)

	nulls=[G_prime.vertex_index[v] for v in G_prime.vertices() if has_neighbor[v]==0]
	if len(nulls)==G.num_vertices():
		return [nulls]
	else:	
		conventions=gen_some_T_conventions(H,T,k,[])
		complements=[node_complements(H,a) for a in conventions]
		cuts = conventions + complements + [nulls]
		return get_intersections(G,cuts)






"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 8: Inference about q
Notes: This section defines the functions used to perform inference about the parameter q. 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""




#Find the share of node v's neighbors adopting the behavior 
def compute_node_q(G,v,ons):
	if G.vertex(v).out_degree()==0:
		return None 

	else: 
		N_v=[G.vertex_index[u] for u in G.vertex(v).out_neighbours()]
		q=len([a for a in N_v if a in ons])/len(N_v)

		return q 




def find_nodewise_q(G,on_nodes,from_prop=False):
	if from_prop is not False:
		ons=[G.vertex_index[v] for v in G.vertices() if from_prop[v]==1]

	else:
		ons=on_nodes

	on_q_list=[]
	off_q_list=[]
	for v in G.vertices():
		q=compute_node_q(G,G.vertex_index[v],ons)
		if q is not None:
			if G.vertex_index[v] in ons:
				on_q_list.append(q)
			else:
				off_q_list.append(q)

	return on_q_list,off_q_list




#For each node in a subset, find the share of its neighbors adopting the behavior
def find_nodewise_q_of_subset(G,on_nodes,subset_prop,subset_prop_val,from_prop=False):
	if from_prop is not False:
		ons=[G.vertex_index[v] for v in G.vertices() if from_prop[v]==1]

	else:
		ons=on_nodes

	on_q_list=[]
	off_q_list=[]
	for v in [v for v in G.vertices() if subset_prop[v]==subset_prop_val]:

		q=compute_node_q(G,G.vertex_index[v],ons)
		if q is not None:
			if G.vertex_index[v] in ons:
				on_q_list.append(q)
			else:
				off_q_list.append(q)

	return on_q_list,off_q_list






#Compare the on-node and off-node histograms
def compare_histograms_q(on_shares,off_shares,num_bins=20,xline=None,height=25,xcolor='#FF00FF',ycolor='#2ECCFA',outfile=None):

	bins = np.linspace(0, 1, num_bins)
	ax = plt.subplot(111)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='on',      # ticks along the bottom edge are off
		top='off',         # ticks along the top edge are off
		labelbottom='on') # labels along the bottom edge are off


	plt.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		left='on',      # ticks along the bottom edge are off
		right='off',         # ticks along the top edge are off
		labelleft='on') # labels along the bottom edge are off


	plt.hist(on_shares, bins, alpha=0.5, label='On Nodes',color=xcolor)
	plt.hist(off_shares, bins, alpha=0.5, label='Off Nodes',color=ycolor)



	if xline is not None:
		plt.axvline(x=xline,linewidth=2,color='#000000',linestyle='--',alpha=0.5)
		plt.text(xline+0.04,height,'$\hat{q}$',color='#000000',alpha=0.5,size=20)

	if outfile is None:
		plt.show()

	else:
		plt.savefig(outfile, transparent=True)


def split_quality_q_of_subset(G,ons,q,subset_prop,subset_prop_val):

	x,y=find_nodewise_q_of_subset(G,ons,subset_prop,subset_prop_val)
	count_1=len([a for a in x if a < q])
	count_2=len([a for a in y if a >= q])
	denom=len(x) + len(y)
	return (1/denom)*(count_1 + count_2)



def split_quality_q(G,ons,q):
	x,y=find_nodewise_q(G,ons)
	count_1=len([a for a in x if a < q])
	count_2=len([a for a in y if a >= q])
	denom=len([v for v in G.vertices() if v.out_degree()>0 or v.in_degree()>0])
	return (1/denom)*(count_1 + count_2)


def find_q(G,ons,include_objval=False,return_interval=False):

	def split_quality_given_G(q):
		return split_quality_q(G,ons,q)

	return minimize_by_line_search(split_quality_given_G,0,1,0.001,include_objval=include_objval,printout=False,return_interval=return_interval)




def find_q_of_subset(G,ons,subset_prop,subset_prop_val,include_objval=False,return_interval=False):

	def split_quality_given_G(q):
		return split_quality_q_of_subset(G,ons,q,subset_prop,subset_prop_val)


	return minimize_by_line_search(split_quality_given_G,0,1,0.001,include_objval=include_objval,printout=False,return_interval=return_interval)





def find_q_from_prop(G,prop,include_objval=False,return_interval=False):
	ons=[G.vertex_index[v] for v in G.vertices() if prop[v]==1]
	if include_objval is False:
		return find_q(G,ons,return_interval=return_interval)
	else:
		return find_q(G,ons,include_objval=True,return_interval=return_interval)



def find_q_of_subset_from_prop(G,prop,subset_prop,subset_prop_val,include_objval=False,return_interval=False):
	ons=[G.vertex_index[v] for v in G.vertices() if prop[v]==1]
	if include_objval is False:
		return find_q_of_subset(G,ons,subset_prop,subset_prop_val,return_interval=return_interval)
	else:
		return find_q_of_subset(G,ons,subset_prop,subset_prop_val,include_objval=True,return_interval=return_interval)




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 9: Inference about t
Notes: This section defines the functions used to perform inference about the parameter t. 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


def compute_node_T(G,v,on):
	if G.vertex(v).out_degree()==0:
		return None 
	else: 
		N_v=[G.vertex_index[u] for u in G.vertex(v).out_neighbours()]
		T=len([a for a in N_v if a in on])

		return T 



def find_nodewise_T(G,ons,from_prop=False):
	if from_prop is False:
		onset=ons
	else:
		onset=[G.vertex_index[v] for v in G.vertices() if from_prop[v]==1]

	on_T_list=[]
	off_T_list=[]
	for v in G.vertices():
		T=compute_node_T(G,G.vertex_index[v],onset)
		if T is not None:
			if v in onset:
				on_T_list.append(T)
			else:
				off_T_list.append(T)

	return on_T_list,off_T_list




def find_nodewise_T_of_subset(G,ons,subset_prop,subset_prop_val,from_prop=False):
	if from_prop is False:
		onset=ons
	else:
		onset=[G.vertex_index[v] for v in G.vertices() if from_prop[v]==1]

	on_T_list=[]
	off_T_list=[]
	for v in  [v for v in G.vertices() if subset_prop[v]==val]:
		T=compute_node_T(G,G.vertex_index[v],onset)
		if T is not None:
			if v in onset:
				on_T_list.append(T)
			else:
				off_T_list.append(T)

	return on_T_list,off_T_list



def split_quality_T(G,ons,T):
	offs=node_complements(G,ons)

	x,y=find_nodewise_T(G,ons)
	count_1=len([a for a in x if a < T])
	count_2=len([a for a in y if a >= T])
	denom=len([v for v in G.vertices() if v.out_degree()>0 or v.in_degree()>0])
	return (1/denom)*(count_1 + count_2)






def find_T(G,ons,include_objval=False):
	offs=node_complements(G,ons)

	consider=list(range(0,max([v.out_degree() for v in G.vertices()])))
	T=0
	working_quality=split_quality_T(G,ons,T) 

	while len(consider)>0:
		newT=consider.pop(0)
		test_quality=split_quality_T(G,ons,newT)
		if test_quality<working_quality:
			T=newT
			working_quality=test_quality

	if include_objval is False:
		return T 
	else:
		return T,working_quality


def find_T_of_subset(G,ons,subset_prop,subset_prop_val,include_objval=False):
	offs=node_complements(G,ons)

	consider=list(range(0,max([v.out_degree() for v in G.vertices()])))
	T=0
	working_quality=split_quality_T(G,ons,T) 

	while len(consider)>0:
		newT=consider.pop(0)
		test_quality=split_quality_T(G,ons,newT)
		if test_quality<working_quality:
			T=newT
			working_quality=test_quality

	if include_objval is False:
		return T 
	else:
		return T,working_quality




def find_T_from_prop(G,prop,include_objval=False):
	ons=[G.vertex_index[v] for v in G.vertices() if prop[v]==1]
	if include_objval is False:
		return find_T(G,ons)
	else:
		return find_T(G,ons,include_objval=True)


def compare_histograms_T(on_shares,off_shares,num_bins=20,xline=None,height=12,xcolor='#FF00FF',ycolor='#2ECCFA',rightshift=0.1,outfile=None,maxT=None):


	minT=min(min(on_shares),min(off_shares)) - 1/2

	if maxT is None:
		maxT=max(max(on_shares),max(off_shares)) + 1/2

	bins = np.arange(minT,maxT,1)
	ax = plt.subplot(111)

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='on',      # ticks along the bottom edge are off
		top='off',         # ticks along the top edge are off
		labelbottom='on') # labels along the bottom edge are off


	plt.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		left='on',      # ticks along the bottom edge are off
		right='off',         # ticks along the top edge are off
		labelleft='on') # labels along the bottom edge are off


	plt.hist(on_shares, bins, alpha=0.5, label='On Nodes',color=xcolor)
	plt.hist(off_shares, bins, alpha=0.5, label='Off Nodes',color=ycolor)



	if xline is not None:
		plt.axvline(x=xline,linewidth=2,color='#000000',linestyle='--',alpha=0.8)
		plt.text(xline+rightshift,height,'$\hat{t}$',color='#000000',alpha=0.8,size=20)

	if outfile is None:
		plt.show()

	else:
		plt.savefig(outfile, transparent=True)





print('Behavioral Communitites functions loaded.')





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 10: Example code
Notes: This section gives an example showing how to find the atomic structure of a random graph
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#Generate an Erdos-Renyi random graph with 30 nodes and density 0.1
G=gen_Erdos_Renyi(30,0.1)


#Find the q=0.5 atomic structure of G using our approximation algorithm 
atoms1=gen_q_range_atoms_approx(G,0.5,3)

#Draw the q=0.5 atoms 
draw_partition(G,atoms2)

#Find the Q=[1/3,1/2]-range atoms 
atoms2=gen_q_range_atoms_approx(G,[1/3,1/2],3)

#Draw the Q=[1/3,1/2]-range atoms 
draw_partition(G,atoms2)


#Find the t=2 atoms 
atoms3=gen_T_atoms_approx(G,2,3)

#Print out the figures showing the atomic stuctures for the atoms computed above 
block_list=[atoms_one_fourth,atoms_one_fourth_robust,atoms_one_half_robust,atoms_one_half_pm_one_tenth,atoms_one_half_pm_one_fifth,atoms_one_third_robust,atoms_one_third_pm_one_tenth,atoms_one_third_pm_one_fifth]

outfile_list=['26KR_atoms_one_fourth.png','26KR_atoms_one_fourth_robust.png','26KR_atoms_one_half_robust.png','26KR_atoms_one_half_pm_one_tenth.png','26KR_atoms_one_half_pm_one_fifth.png','26KR_atoms_one_third_robust.png','26KR_atoms_one_third_pm_one_tenth.png','26KR_atoms_one_third_pm_one_fifth.png']

draw_comparisons_cc(G,block_list,outfile_list,shape=G.vp.caste,vertex_halo=False)





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 12 - AddHealth data analysis 

Note: the AddHealth data is not public and cannot be shared, and so the datafiles used in the proceeding code are not included  
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

dataset=pd.read_csv('AddHealth Data/network_behavior.csv',dtype={'aid':str,'mf1aid':str,'mf2aid':str,'mf3aid':str,'mf4aid':str,'mf5aid':str,'ff1aid':str,'ff2aid':str,'ff3aid':str,'ff4aid':str,'ff5aid':str})


#Look at sister schools 91 and 191
dataset=dataset.loc[dataset['sschlcde'].isin([91,191])]

#dataset=dataset.loc[dataset['sschlcde'].isin([44,144])]
#dataset=dataset.loc[dataset['sschlcde'].isin([31])]


#Drop observations missing important data 

dataset=dataset.loc[dataset['aid'].notnull() & dataset['s59a'].notnull()]
#dataset=dataset.loc[dataset['aid'].notnull() & dataset['s49'].notnull()]

dataset.sort_values('aid')


for var in ['s2','s3','s4','s6a','s6b','s6c','s6d','s1','s5','s49','s59a']:

  for i in dataset.index.values:
      if isnan(dataset[var][i]):
        dataset.set_value(i,var,999)
  

  dataset[var]=dataset[var].astype(int)





dataset['race']=str('what')

#Code the race as 1=black, 2=hispanic, 3=Asian, 4=White, 5=Other 
for i in dataset.index.values:
  if dataset['s4'][i]==1:
    dataset.set_value(i,'race','hispanic')
  else: 
    if dataset['s6b'][i]==1:
      dataset.set_value(i,'race','black')
    elif dataset['s6c'][i]==1:
      dataset.set_value(i,'race','Asian')
    elif dataset['s6a'][i]==1:
      dataset.set_value(i,'race','white')
    else:
      dataset.set_value(i,'race','other')



#dataset=dataset.loc[dataset['sschlcde'].isin([91,191])]


#Convert the dataset to a network 

dataset=dataset.reset_index(drop=True)

node_labels=pd.Series([int(a) for a in range(dataset.shape[0])])
dataset['node_label']=node_labels


node_dict=dict(zip(list(dataset['aid']),list(dataset['node_label'])))
smoke_dict=dict(zip(list(dataset['node_label']),[int(a) not in [0,99] for a in list(dataset['s59a'])]))
sex_dict=dict(zip(list(dataset['node_label']),[int(a) >1 for a in list(dataset['s2'])]))
grade_dict=dict(zip(list(dataset['node_label']),[int(a) for a in list(dataset['s3'])]))
race_dict=dict(zip(list(dataset['node_label']),[a for a in list(dataset['race'])]))
drink_dict=dict(zip(list(dataset['node_label']),[int(a) not in [0,99] for a in list(dataset['s49'])]))

friendlistvars=['mf'+str(i)+'aid' for i in range(1,6)] + ['ff'+str(i)+'aid' for i in range(1,6)]

G=Graph()
G.add_vertex(dataset.shape[0])


for i in range(dataset.shape[0]):
  for var in friendlistvars:
    entry=dataset[var][i]
    if entry in list(node_dict.keys()):
      neighbor_node=node_dict[entry]
      G.add_edge(i,neighbor_node)


#G=convert_to_undirected(G,reciprocal_only=False)
G=convert_to_undirected(G,reciprocal_only=True)
G=get_big_component(G)




draw_with_colors(G,component)

gen_node_property(G,smoke_dict,'smokes')
gen_node_property(G,sex_dict,'isgirl')
gen_node_property(G,grade_dict,'grade')
gen_node_property(G,race_dict,'race',type='string')
gen_node_property(G,drink_dict,'drinks')


gen_node_property(G_reciprocal,smoke_dict,'smokes')
gen_node_property(G_reciprocal,sex_dict,'isgirl')
gen_node_property(G_reciprocal,grade_dict,'grade')
gen_node_property(G_reciprocal,race_dict,'race',type='string')
gen_node_property(G_reciprocal,drink_dict,'drinks')



error_rates=[]
atom_list=[]

for q in [0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
  atoms=gen_q_atoms_approx_shortcut_multicore(G,q,3)
  error_rates.append(partition_error_rate(G,G.vp.smokes,atoms))
  atom_list.append(atoms)




dwcc(G,G.vp.smokes,vertex_size=10,outfile='smokes_allboys.png')

q_est,q_val=find_q_from_prop(G,G.vp.smokes,include_objval=True,return_interval=True)
T_est,T_val=find_T_from_prop(G,G.vp.smokes,include_objval=True)

q_point_est=mean(q_est)

ons,offs=find_nodewise_q(G,[],from_prop=G.vp.smokes)
compare_histograms_q(ons,offs,num_bins=12,xline=q_point_est,height=120)
compare_histograms_q(ons,offs,num_bins=15,xline=q_point_est,height=120,outfile='addhealth_q_split_allboys.png')

#compare_histograms_q(ons,offs,xline=q_est,height=120,outfile='addhealth_q_split.png')


compare_histograms_q(ons,offs,xline=q_est,height=120)


ons,offs=find_nodewise_T(G,[],from_prop=G.vp.smokes)
#compare_histograms_T(ons,offs,xline=T_est,height=140,rightshift=0.9,outfile='addhealth_t_split.png')
compare_histograms_T(ons,offs,xline=T_est,height=140,num_bins=15,rightshift=0.9,outfile='addhealth_t_split_allboys.png',maxT=12)



q_est_girls,q_val_girls=find_q_of_subset_from_prop(G,G.vp.smokes,subset_prop=G.vp.isgirl,subset_prop_val=1,include_objval=True)
q_est_boys,q_val_boys=find_q_of_subset_from_prop(G,G.vp.smokes,subset_prop=G.vp.isgirl,subset_prop_val=0,include_objval=True)


ons,offs=find_nodewise_q_of_subset(G,[],subset_prop=G.vp.isgirl,subset_prop_val=1,from_prop=G.vp.smokes)
compare_histograms_q(ons,offs,xline=q_est_girls,num_bins=10,height=100,outfile='addhealth_q_split_girls.png')

ons,offs=find_nodewise_q_of_subset(G,[],subset_prop=G.vp.isgirl,subset_prop_val=0,from_prop=G.vp.smokes)
compare_histograms_q(ons,offs,xline=q_est_boys,num_bins=10,height=100,outfile='addhealth_q_split_boys.png')



grade_results=[]
for grade in [9,10,11,12]:
  grade_results.append([grade,find_q_of_subset_from_prop(G,G.vp.smokes,subset_prop=G.vp.grade,subset_prop_val=grade,include_objval=True,return_interval=True)])


counter=0
for grade in [9,10,11,12]:

  ons,offs=find_nodewise_q_of_subset(G,[],subset_prop=G.vp.grade,subset_prop_val=grade,from_prop=G.vp.smokes)
  compare_histograms_q(ons,offs,xline=mean(grade_results[counter][1][0]),num_bins=12,height=30,outfile='addhealth_q_split_grade_'+str(grade)+'_allboys.png')
  plt.close() 
  counter=counter+1





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 13 - Zacharys Karate Club 

Note: This section includes the code to reproduce figures 25-27 in the paper  
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

G=load_graph('karate.gt')
loyal=[G.vertex_index[v] for v in G.vertices() if loyalty[v] ==0]

block_list=[[],[loyal]]
outfile_list=['karate_blank.png','karate_club.png']

#Draw Zachary's Karate club with node color indicating whether a node joined the splinter group (Figure 25)
draw_comparisons_cc(G,color_list=['#2ECCFA','#FF00FF'],block_list=block_list,vertex_size=30, outfile_list=outfile_list,output_size=1500)

disloyal=[G.vertex_index[v] for v in G.vertices() if loyalty[v]==1]

x,y=find_nodewise_q(G,disloyal)


#Compare the on-node and off-node on-neighbor shares histograms (Figure 26)
compare_histograms(x,y,num_bins=10,height=12,outfile='karate_hist_a.png')
time.sleep(3)
plt.close()
time.sleep(3)
compare_histograms(x,y,num_bins=10,height=12,xline=0.43,outfile='karate_hist_b.png')
time.sleep(3)
plt.close()





#Compare how well different choices of q split the on-node and off-node shares (Figure 27) 
split_qualitites=[]
for q in np.arange(0.1,0.9,0.05):
    split_qualitites.append(split_quality(G,disloyal,q))

line_plot([np.arange(0.1,0.9,0.05)],[split_qualitites])
scatter_plot(np.arange(0.1,0.9,0.05),split_qualitites,marker_size=100)



scatter_plot(np.arange(0.1,0.9,0.05),split_qualitites,xlabel='potential q',ylabel='T(q)',marker_size=100,outfile='exploring_karate_splits.png')




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 14 - Seeding Simulations 

Note: This section reproduces the seeding simulations from section 6 of the paper   
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#Load collection of 10 village graphs

os.chdir('Accompanying Code and Data (Third Submission)')
village_HH_kr_graphs=[]

for i in (set(range(1,78)) - {13,22}):
    file_to_read='Indian Villages HH/Networks/adj_keroricego_HH_vilno_' + str(i)  + '.csv'
    reader = csv.reader(open(file_to_read, "r"), delimiter=",")
    x = list(reader)
    result = np.array(x).astype("int")
    G=adj_to_graph(result)
    village_HH_kr_graphs.append(G)


G_set=village_HH_kr_graphs[6, 12, 14, 27, 28, 44, 50, 53, 65, 74]



#Loop over q in a range and compare the seeding results for random, optimal, and atom-based heuristic seeding for Indian village number 77

q_list=[a/100 for a in range(20,81,10)]
num_seeds=5 

rand_spread=[]
heur_spread=[]
opt_spread=[]
Lou_spread=[]


for q in q_list:

    for G in G_set:

    rand_temp=[]
    for k in range(0,1000):
        rand_temp.append(random_seed_result(G,q,num_seeds))
    rand_spread.append(mean(rand_temp))
	opt_spread.append(optimal_seeding_by_brute_force(G,q,num_seeds)[0])
	heur_spread.append(atom_based_heuristic(G,q,num_seeds)[0])
	Lou_spread.append(Louvain_based_heuristic(G,q,num_seeds))


rand_spread=[a/155 for a in rand_spread]
heur_spread=[a/155 for a in heur_spread]
opt_spread=[a/155 for a in opt_spread]




"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Section 2: CalTech Data Analysis (Analysis of Link Persistence in Caltech Data)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#Load the Caltech networks 
dataset = pd.read_stata('caltechData/DyadWave2013Dataset_noR_newObsMarked.dta')
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100)


dataset15=dataset[dataset.wave=='F15']
node_ids15=dataset15.from_uid.unique()

dataset14=dataset[dataset.wave=='F14']
node_ids14=dataset.from_uid.unique()

node_ids=list(set(node_ids14).intersection(set(node_ids15)))

dataset14=dataset14[dataset14.to_uid.isin(node_ids)]
dataset14=dataset14[dataset14.from_uid.isin(node_ids)]

dataset15=dataset15[dataset15.to_uid.isin(node_ids)]
dataset15=dataset15[dataset15.from_uid.isin(node_ids)]

G=Graph(directed=True)

G.add_vertex(len(node_ids))

id_dict={}

count=0
for v in G.vertices():
    id_dict[v]=node_ids[count]
    count=count+1

gen_node_property(G,id_dict,'id')

inv_map = dict(zip(id_dict.values(), id_dict.keys()))

G14=G.copy()
G15=G.copy()


links14=[]

dataset14=dataset14.reset_index(drop=True)

for i in range(dataset14.shape[0]):
    if dataset14.loc[i].FriendDirected==1:
        links14.append([inv_map[dataset14.loc[i].from_uid],inv_map[dataset14.loc[i].to_uid]])

G14.add_edge_list(links14)
G14=convert_to_undirected(G14,reciprocal_only=False)

links15=[]

dataset15=dataset15.reset_index(drop=True)

for i in range(dataset15.shape[0]):
    if dataset15.loc[i].FriendDirected==1:
        links15.append([inv_map[dataset15.loc[i].from_uid],inv_map[dataset15.loc[i].to_uid]])

G15.add_edge_list(links15)
G15=convert_to_undirected(G15,reciprocal_only=False)





compare_graphs(G14,G15,space=15)

atoms14=gen_q_range_atoms_approx_multicore(G14,[1/3- epsilon, 1/3 + epsilon],3)
atoms14alt=gen_q_atoms_approx_shortcut_multicore(G14,1/3,3)



make_atoms_vp(G14,atoms14)

within_atom_edges=[]
across_atom_edges=[]
persistent_edges=[]



for e in G14.edges():
    source_index=G14.vertex_index[e.source()]
    target_index=G14.vertex_index[e.target()]

    if G14.vp.atom_number[source_index]==G14.vp.atom_number[target_index]:
        within_atom_edges.append(e)
    else:
        across_atom_edges.append(e)

    if G15.edge(source_index,target_index) is not None:
        persistent_edges.append(e)


print(len(set(persistent_edges).intersection(set(within_atom_edges)))/len(within_atom_edges))
print(len(set(persistent_edges).intersection(set(across_atom_edges)))/len(across_atom_edges))
