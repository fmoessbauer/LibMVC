#include "tsewf.h"

const int try_step = 600000;

inline unsigned int rdrand32 ()
{
	register unsigned int rand;
	asm volatile ("rdrand %0"
		: "=r" (rand)
	);
	return rand;
}

void forget_edge_weights()
{
	int v,e;
	int new_total_weight=0;
	
	for(v=1; v<=v_num; v++)
		dscore[v]=0;

	//scale_ave=ave_weight*q_scale;
	for (e = 0; e<e_num; e++)
	{
		edge_weight[e] = edge_weight[e]*p_scale;
			
		new_total_weight+=edge_weight[e];
		
		//update dscore
		if (v_in_c[edge[e].v1]+v_in_c[edge[e].v2]==0){
			dscore[edge[e].v1]+=edge_weight[e];
			dscore[edge[e].v2]+=edge_weight[e];
			}
		else if(v_in_c[edge[e].v1]+v_in_c[edge[e].v2]==1){
			if(v_in_c[edge[e].v1]==1)dscore[edge[e].v1]-=edge_weight[e];
			else  dscore[edge[e].v2]-=edge_weight[e];
		}
	}
	ave_weight=new_total_weight/e_num;

}


void update_edge_weight()
{
	int i,e;
	for(i=0; i<uncov_stack_fill_pointer; ++i)
	{
		e = uncov_stack[i];
		
		edge_weight[e]+= 1;
		dscore[edge[e].v1] += 1;
		dscore[edge[e].v2] += 1;
	}
	
	
	delta_total_weight += uncov_stack_fill_pointer;
	if(delta_total_weight >= e_num)
	{
		ave_weight += 1;
		delta_total_weight -= e_num;
	}
	
	//smooth weights
	if(ave_weight >= threshold)
	{
		forget_edge_weights();
	}
	
}


void cover_LS()
{
	int		best_cov_v;		//the vertex of the highest dscore in C
	int		best_add_v;
	int		e,v1,v2;
	int		i,v;

	step = 1;

	while(1)// wihin cutoff_time
	//while(step<=max_steps)
	{
		/* ### if there is no uncovered edge ### */
		if (uncov_stack_fill_pointer == 0)
		{
			update_best_sol();			// C* := C
			
			printf("Better MVC found.\tSize: %d\tTime: %.3f\n", v_num-best_c_size,best_comp_time);
			
			if (c_size==optimal_size)
				return;
				
			update_target_size();		// remove a vertex with the highest dscore from C;
			
			continue;
		}
		
		/* monitor status */
		
		if(step % try_step==0)
		{
			finish = clock();
			double elap_time = ((double)finish - start_time)/CLOCKS_PER_SEC;
			//printf("Best MVC size: %d\tTime: %.3f\n",best_c_size,elap_time);
			if(elap_time >= cutoff_time)return;
		}
		
		/* choose a vertex u in C with the highest dscore,
			breaking ties in favor of the oldest one; */
		best_cov_v = remove_cand[0];
		for (i=1; i<remove_cand_size; ++i)
		{
			v = remove_cand[i];
			if(v==tabu_remove)
				continue;
			if( dscore[v] < dscore[best_cov_v])
				continue;
			else if( dscore[v] > dscore[best_cov_v] )
				best_cov_v = v;
			else if (time_stamp[v] < time_stamp[best_cov_v])
				best_cov_v = v;
		}
	
		/* C := C\{u}, confChange(u) := 0 and confChange(z) := 1 for each z in N(u); */
		remove(best_cov_v);
		
		
		/* choose an uncovered edge e randomly; */
		//e = uncov_stack[rand()%uncov_stack_fill_pointer];
		e = uncov_stack[rdrand32()%uncov_stack_fill_pointer];
		
		/* choose a vertex v in e such that confChange(v) = 1 with higher dscore, 
			breaking ties in favor of the older one; */
		v1 = edge[e].v1;
		v2 = edge[e].v2;
		
		if(conf_change[v1]==0 ) best_add_v=v2;
		else if(conf_change[v2]==0) best_add_v=v1;
		else
		{
			if(dscore[v1]>dscore[v2] || (dscore[v1]==dscore[v2] && time_stamp[v1]<time_stamp[v2]) )
				best_add_v=v1;
			else
				best_add_v=v2;
		}
		
		
		/* C := C plus {v}, confChange(z) := 1 for each z in N(v); */
		add(best_add_v);
		
		// update remove candidate index. speed optimization.
		int index = index_in_remove_cand[best_cov_v];
		index_in_remove_cand[best_cov_v] = 0;
		
		remove_cand[index] = best_add_v;
		index_in_remove_cand[best_add_v] = index;
		
		// update timestamp used to determine the age of each virtex
		time_stamp[best_add_v] = step;
		time_stamp[best_cov_v] = step;
		
		tabu_remove = best_add_v;
		
		/* w(e) := w(e) + 1 for each uncovered edge e; */
		/* if w >= y then w(e) := [p*w(e)] for each edge e; */
		update_edge_weight();
		
		step++;
	}
}

int main(int argc, char* argv[])
{
	int seed;
	char title[1024];
	
	if(build_instance(argv[1])!=1){
		cout<<"can't open instance file"<<endl;
		return -1;
	}
	
		sscanf(argv[2],"%d",&optimal_size);//if you want to stop the algorithm only cutoff time is reached, set optimal_size to 0.
		sscanf(argv[3],"%d",&seed);
		sscanf(argv[4],"%d",&cutoff_time);

	sprintf(title,"@title file: %s        optimal_size: %d        seed: %d",argv[1],v_num-optimal_size,seed);
	system(title);
	/*
	if(build_instance("frb45-21-1.mis")!=1){
		cout<<"can't open instance file"<<endl;
		return -1;
	}
		optimal_size = 900;
		seed = 1100;
		cutoff_time = 1000;
		
	*/
	
		threshold = (int)(0.5*v_num); 
	
	
		srand(seed);

		//cout<<seed<<' ';
		cout<<"c Improved NuMVC Local Search Solver"<<endl;

    	init_sol();
    	
    	start = clock();
		start_time = start;
		
	if(c_size + uncov_stack_fill_pointer > optimal_size ) 
	{
		cout<<"c Start local search..."<<endl;
		cover_LS();
	}
		
		//check solution
		if(check_solution()==1)
		{
			cout<<"c Best found vertex cover size = "<<best_c_size<<endl;
			print_solution();
			cout<<"c searchSteps = "<<best_step<<endl;
			cout<<"c solveTime = "<<best_comp_time<<endl;
			cout<<"c performance = "<<(best_step/best_comp_time)/1000000<< "MT/s"<<endl;
			
			//cout<<best_c_size<<' '<<best_comp_time<<' '<<best_step<<endl;
		}
		else
		{
			cout<<"the solution is wrong."<<endl;
			//print_solution();
		}
	
		free_memory();

	return 0;
}
