
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>

/* global parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,105,110,115,120, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310};
int NUM_OF_RUNS = 1;
int MAX_TIME = 300;  //max amount of time permited (in sec)
int num_of_problems;

/* parameters for evlutionary algorithms*/
static int POP_SIZE = 10;   //please modify these parameters according to your problem
int MAX_NUM_OF_GEN = 10000; //max number of generations
float MUTATION_RATE = 0.003;
float REPLACE_RATE = 0.8;

/* declare parameters for variable neighbourhood search here*/
int K= 3; // k-opt is used

//return a random number between 0 and 1
int rand_01()
{

    int number;
    number = rand() % 2;
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}

//struct of item
struct item_struct{

    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
    double ratio;
    int indx;

};

//struct of problem
struct problem_struct{

    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
    double best_obj; //best known objective

};

//free the problem
void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

//init the problem
void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}

//load problem
struct problem_struct** load_problems(char* data_file)
{
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);
 
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n);
        fscanf (pfile, "%d", &dim); fscanf (pfile, "%d", &obj_opt);
        
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;
            my_problems[k]->items[j].indx=j;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
            //printf("item[j].p=%d\n",my_problems[k]->items[j].p);
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
                //printf("my_problems[%i]->items[%i].size[%i]=%d\n",k,j,i,my_problems[k]->items[j].size[i]);
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
            //printf("capacities[i]=%d\n",my_problems[k]->capacities[i] );
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

//struct of solution
struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;
    int feasibility; //indicate the feasiblity of the solution
    int* x; // solution encoding vector
    int* cap_left; //capacity left in all dimensions
};

//init the solution
struct solution_struct init_solution(struct problem_struct* my_prob)
{
    int genenum = my_prob->n;
    int dimnum = my_prob->dim;
    struct solution_struct* new_sln = malloc(sizeof(struct solution_struct));
    new_sln->prob  = my_prob;
    new_sln->objective = 0;
    new_sln->feasibility = 1;
    new_sln->x = malloc(sizeof(int)*genenum);
    new_sln->cap_left = malloc(sizeof(int)*dimnum);
    return *new_sln;
}

//free the solution
void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
    }
}

//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

//evaluate the solution
void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution
    sln->objective =0; sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;
    
    for(int i=0; i< items_p->dim; i++)
    {
        sln->cap_left[i]=sln->prob->capacities[i];
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //exceeding capacity
                return;
            }
        }
    }
    if(sln->feasibility>0)
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}


//output a given solution to a file
void output_solution(struct solution_struct* sln, char* out_file)
{
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        fprintf(pfile, "%i\n", (int)sln->objective);
        for(int i=0; i<sln->prob->n; i++)
        {
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}

//randomly insert the gene
void insert_gene(struct solution_struct* populations,int gene_size){
    
    for(int i = 0; i < POP_SIZE; i++) {
	for(int j = 0; j < gene_size; j++){
	    int rand_num = rand_int(0,gene_size-1);
	    populations[i].x[rand_num] = 1;
	    evaluate_solution(&populations[i]);
	    if(populations[i].feasibility <= 0){
		populations[i].x[rand_num] = 0;
		evaluate_solution(&populations[i]);
		break;
	    }
	}
    }

}

//init_the population 
struct solution_struct* init_populations(struct problem_struct* prob){

    int genesize = prob->n;
    struct solution_struct* population= malloc(sizeof(struct solution_struct)*POP_SIZE);
    for (int i = 0; i < POP_SIZE; i++)
    {
	//init each unit
        population[i] = init_solution(prob);
	evaluate_solution(&population[i]);
    }
    insert_gene(population,genesize);
    
    return population;

}

// Random select the parent by using corona
struct solution_struct select_parent(struct solution_struct* population, struct problem_struct* prob){
    int genenum = population[0].prob->n;
    int dimnum = population[0].prob->dim;
    
    //calculate the whole fittness of the population
    float total_ob = 0;
    for (int i = 0; i < POP_SIZE; i++){
	total_ob += population[i].objective;
    }
    
    //create the corona
    float *corona = malloc (sizeof(float)*POP_SIZE*2);
    corona[0] = 0;
    float total_pbb;
    for (int i = 1; i <= POP_SIZE; i++)
    {
	total_pbb += population[i-1].objective/total_ob;
        corona[i] = total_pbb;
    }
    
    //select the parent
    int index = 0;
    float rand_pbb = 1.0 * rand()/RAND_MAX;
    for (int i = 1; i <= POP_SIZE; i++)
    {
	if (corona[i] >= rand_pbb){
	    index = i-1;
	    break;
	}
    }

    //return the parent
    struct solution_struct parent = init_solution(prob);
    copy_solution(&parent, &population[index]);
    free(corona);
    return parent;

}

//set up the mating pool
struct solution_struct* set_mating_pool(struct solution_struct* population,struct problem_struct* prob){

    // init the mating pool
    struct solution_struct* mating_pool;
    mating_pool= malloc(sizeof(struct solution_struct)*POP_SIZE);

    int genenum = population[0].prob->n;
    int dimnum = population[0].prob->dim;
    struct solution_struct parent_temp = init_solution(prob);
    

    for (int i = 0; i < POP_SIZE; i++)
    {
        
        mating_pool[i] = init_solution(prob);
	parent_temp = select_parent(population,prob);
	copy_solution(&mating_pool[i], &parent_temp);

    }

    free_solution(&parent_temp);
    return mating_pool;

}

//min to max
int cmpfunc(const void* a, const void* b){
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->ratio>item2->ratio) return 1;
    if(item1->ratio<item2->ratio) return -1;
    return 0;
}

//max to min
int cmpfunc1(const void* a, const void* b){
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->ratio>item2->ratio) return -1;
    if(item1->ratio<item2->ratio) return 1;
    return 0;
}

//back to init
int cmpfunc2 (const void * a, const void * b) {
        const struct item_struct* item1 = a;
        const struct item_struct* item2 = b;
        if(item1->indx>item2->indx) return 1;
        if(item1->indx<item2->indx) return -1;
        return 0;
}

//min to max
int cmpfunc_sln (const void * a, const void * b) {
    const struct solution_struct* sln1 = a;
    const struct solution_struct* sln2 = b;
    if(sln1->objective > sln2 ->objective) return 1;
    if(sln1->objective < sln2 ->objective) return -1;
    return 0;
}

//max to min
int cmpfunc_sln1 (const void * a, const void * b) {
    const struct solution_struct* sln1 = a;
    const struct solution_struct* sln2 = b;
    if(sln1->objective > sln2 ->objective) return -1;
    if(sln1->objective < sln2 ->objective) return 1;
    return 0;
}

// using greedy to delete the min ratio item until the chile is feasibile
struct solution_struct repair_solution (struct solution_struct child, struct problem_struct* prob){

    int genenum = prob->n;
    int dimnum = prob->dim;

    //calculate the item ratio
    for(int i=0; i<genenum;i++){
        double avg_size=0;
        struct item_struct* item_i = &prob->items[i];
        for(int d=0; d< dimnum; d++){
            avg_size += (double)item_i->size[d]/prob->capacities[d];
        }
        item_i->ratio = item_i->p/avg_size;
    }
    
    //sort from min to max
    qsort(prob->items, genenum, sizeof(struct item_struct), cmpfunc);

    //delete the min ratio of iten until fit the size
    for(int i = 0; i<genenum; i++){

	if(child.x[prob->items[i].indx] == 1){
	    child.x[prob->items[i].indx] = 0;
	    evaluate_solution(&child);
	}
	if (child.feasibility == 1){
	    break;
	}
	
    }

    //sort to normal
    qsort(prob->items, genenum, sizeof(struct item_struct), cmpfunc2);
    return child;
}

struct solution_struct mate (struct solution_struct* mating_pool, struct problem_struct* prob){
    int genenum = prob->n;
    int dimnum = prob->dim;

    //init the solution
    struct solution_struct parent_1 = init_solution(prob);
    struct solution_struct parent_2 = init_solution(prob);
    struct solution_struct child = init_solution(prob);

    //random select the parents in mating pool
    int rand_par1 = rand_int(0,POP_SIZE-1);
    int rand_par2 = rand_int(0,POP_SIZE-1);
    copy_solution(&parent_1,&mating_pool[rand_par1]);
    copy_solution(&parent_2,&mating_pool[rand_par2]);

    //using uniform crossover and mutation
    for(int i = 0; i < genenum; i++){
	float rand_cro = 1.0 * rand()/RAND_MAX;
        float rand_mut = 1.0 * rand()/RAND_MAX;
	if(rand_cro < 0.5){
	    child.x[i] = parent_1.x[i];
	}
	else{
	    child.x[i] = parent_2.x[i];
	}
	// mutation
        if(rand_mut <= MUTATION_RATE){
	    if(child.x[i] == 0){
		child.x[i] = 1;
	    }
	    else{
		child.x[i] = 0;
	    }
	}
    }

    free_solution(&parent_1);
    free_solution(&parent_2);
    evaluate_solution(&child);

    if (child.feasibility != 1){
	    struct solution_struct child_fix = repair_solution(child, prob);
        return child_fix;
    }
    else{
        return child;
    }
}

bool can_swap(struct solution_struct* sln, int out, int in)
{
    for(int d =0; d<sln->prob->dim; d++)
    {
        if(sln->cap_left[d]+sln->prob->items[out].size[d] < sln->prob->items[in].size[d])
            return false;
    }
    return true;
}



bool can_move(int nb_indx, int* move, struct solution_struct* curt_sln ){
    bool ret=true;
    if(nb_indx==1)
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<curt_sln->prob->dim; d++){
            if(curt_sln->cap_left[d] < curt_sln->prob->items[i].size[d])
                return false;
        }
    }
    else if(nb_indx==2){
        ret=can_swap(curt_sln, move[0], move[1]);
    }
    else if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(curt_sln->x[j]>0) {//2-1 swap
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] +
                   curt_sln->prob->items[j].size[d] < curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        else {//1-2 swap
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] <
                   curt_sln->prob->items[j].size[d] + curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        
    }
    else ret=false;
    return ret;
}

bool apply_move(int nb_indx, int* move, struct solution_struct* sln ){
    bool ret=true;
    if(nb_indx==1)
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] -= sln->prob->items[i].size[d];
        }
        sln->objective += sln->prob->items[i].p;
        sln->x[i]=1;
        
        //printf("success\n");
    }
    else if(nb_indx==2){
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] = sln->cap_left[d] + sln->prob->items[move[0]].size[d]-
                sln->prob->items[move[1]].size[d];
        }
        sln->objective += sln->prob->items[move[1]].p-sln->prob->items[move[0]].p;
        sln->x[move[0]]=0; sln->x[move[1]]=1;
    }
    else if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(sln->x[j]>0) {//2-1 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] +
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[k].p-sln->prob->items[i].p-sln->prob->items[j].p;
            sln->x[i]=0; sln->x[j]=0; sln->x[k]=1;
        }
        else {//1-2 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] -
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[j].p+sln->prob->items[k].p-sln->prob->items[i].p;
            sln->x[i]=0; sln->x[j]=1; sln->x[k]=1;
        }
        
    }
    else ret=false;
    return ret;
}

//nb_indx <=3
struct solution_struct* best_descent_vns(int nb_indx, struct solution_struct* curt_sln)
{
    struct solution_struct* best_neighb = malloc(sizeof(struct solution_struct));
    best_neighb->cap_left = malloc(sizeof(int)*curt_sln->prob->dim);
    best_neighb->x = malloc(sizeof(int)*curt_sln->prob->n);
    copy_solution(best_neighb, curt_sln);
    int n=curt_sln->prob->n;
    int curt_move[] ={-1,-1,-1}, best_move []={-1,-1,-1}, delta=0, best_delta=0;  //storing best neighbourhood moves
    switch (nb_indx)
    {
        case 1: //check whether any items can be inserted.
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]>0) continue;
                curt_move[0]=i;
                if(can_move(nb_indx, &curt_move[0], best_neighb)){
                    delta = curt_sln->prob->items[i].p;
                    if(delta> best_delta) {
                        best_delta = delta; best_move[0] = i;
                    }
                }
            }
            if(best_delta>0) {    apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 2:
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]<=0) continue;
                for(int j=0; j<n; j++){
                    if(curt_sln->x[j]==0)
                    {
                        curt_move[0]= i; curt_move[1]= j; curt_move[2]=-1;
                        if(can_move(nb_indx, &curt_move[0], best_neighb)){
                            delta = curt_sln->prob->items[j].p -curt_sln->prob->items[i].p;
                            if(delta > best_delta){
                                best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=-1;
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 3:
            //2-1 swap
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j!=i&&j<n; j++){
                    if(curt_sln->x[j]==0) continue;
                    for(int k=0;k<n;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], best_neighb)){
                                delta = curt_sln->prob->items[k].p -curt_sln->prob->items[i].p-curt_sln->prob->items[j].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                }
                            }
                        }
                    }
                }
            }
            //1-2 swap
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j<n; j++){
                    if(curt_sln->x[j]>0) continue;
                    for(int k=0;k!=j&&k<n;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], curt_sln)){
                                delta = curt_sln->prob->items[k].p +curt_sln->prob->items[j].p-curt_sln->prob->items[i].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                }
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        default:
            printf("Neighbourhood index is out of the bounds, nothing is done!\n");
    }
    return best_neighb;
}


//VNS
void varaible_neighbourhood_search(struct solution_struct* my_sln, struct problem_struct* prob, clock_t time_start){

    int nb_indx =0; //neighbourhood index
    struct solution_struct curt_sln = init_solution(prob);
    struct solution_struct* neighb_s;
    clock_t time_fin;
    double time_spent=0;
    for(int i = 0; i < POP_SIZE; i++){
        copy_solution(&curt_sln,&my_sln[i]);
    
	    int count = 0;
        while(count < 4){
            while(nb_indx<K){
                neighb_s = best_descent_vns(nb_indx+1, &curt_sln); //best solution in neighbourhood nb_indx
                if(neighb_s->objective > curt_sln.objective){
                    copy_solution(&curt_sln, neighb_s);
                    nb_indx=1;
                }
                else nb_indx++;
            }
            count++;
            nb_indx=0;
        }
        copy_solution(&my_sln[i],&curt_sln);
 	    //if time > MAX TIME out
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
        if(time_spent >= MAX_TIME-10){break;}
    }
    
    free_solution(neighb_s);free(neighb_s);
    free_solution(&curt_sln);

}



int MA(struct problem_struct* prob, char* out_file){

    struct solution_struct* my_pop = init_populations(prob);
    struct solution_struct* mating_pool;
    struct solution_struct child = init_solution(prob);
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    
    
    for(int i = 0; i < 500; i++){
        //if time > MAX TIME out
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
        if(time_spent >= MAX_TIME-10){goto end;}
        for(int k = 0; k < 200; k++){
            mating_pool = set_mating_pool(my_pop,prob);
            child = mate(mating_pool,prob);
            qsort(my_pop, POP_SIZE, sizeof(struct solution_struct), cmpfunc_sln);
            //Generation gap
            for(int j=0; j<POP_SIZE*REPLACE_RATE; j++){
                copy_solution(&my_pop[j],&child);
            }

            //if time > MAX TIME out
            time_fin=clock();
            time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
            if(time_spent >= MAX_TIME-10){goto end;}
        }

        qsort(my_pop, POP_SIZE, sizeof(struct solution_struct), cmpfunc_sln1);

        varaible_neighbourhood_search(my_pop, prob, time_start);
    }
	
    end:
    

    //double check the feasibility
    for(int i = 0; i<POP_SIZE; i++){
	evaluate_solution(&my_pop[i]);
    }

    //sort the my_pop
    qsort(my_pop, POP_SIZE, sizeof(struct solution_struct), cmpfunc_sln1);
    output_solution(&my_pop[0],out_file);
    printf("ob: %f\n", my_pop[0]. objective);
    
    for(int k=0; k<POP_SIZE; k++)
    {
        free_solution(&my_pop[k]); //free pop and pool data memory
        free_solution(&mating_pool[k]);
    }
    
    free_solution(&child);
    free(my_pop);
    free(mating_pool);
    return 0;

}

int main(int argc, const char * argv[]) {
    
    printf("Starting the run!\n");
    char data_file[50]={"somefile"}, out_file[50]={}, solution_file[50]={};  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        //printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    }
    
    struct problem_struct** my_problems = load_problems(data_file);

    FILE* pfile = fopen(out_file, "w"); //open a new file
    fprintf(pfile, "%d\n", num_of_problems); fclose(pfile);
    for (int i = 0; i < num_of_problems; i++){
        printf("Problem %d, running\n",i);
        srand(RAND_SEED[i%30]);
	    MA(my_problems[i],out_file);
    }
    

    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }

    free(my_problems); //free problems array
    
    return 0;
}
