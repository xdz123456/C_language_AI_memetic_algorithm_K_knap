**1.**  **Description**

Memetic algorithm is a combination of local and global search algorithm. In my program, I select the genetic algorithm as the global search and the VNS algorithm as the local search. 

 

**1.**   **Initial:**

Firstly, I initialized the population and randomly selected items until any dimension of the answer exceeded the size of the knapsack. Algorithm 1 is the pseudo code of initial Population.

**Algorithm 1: Initial Population****：**

Input: The current problem Prob 

Output: The initialized population P which a set of solution

```
\1.   Define population P with size of POP_SIZE



\2.   Genesize = Prob.n

\3.   for i = 0 to POP_SIZE - 1 do

\4.     for j = 0 to to Genesize - 1 do

\5.         P[i].x[j] = 0

\6.     end for

\7.     for k = 0 to Genesize - 1 do

\8.       generate random index rand_index

\9.       P[i].x[rand_index] = 1

\10.      evaluate_solution(P[i])

\11.      if(P[i].feasibility <= 0)

\12.        P[i].x[rand_index] = 0

\13.        evaluate_solution(P[i])

\14.        break

\15.      end if

\16.     end for

\17.  end for

\18.  return P
```

**2.**   **Selection:**

I select parent generation through Roulette Wheel Selection and establish a mating pool which size is POP_SIZE, Algorithm 2 is the pseudo code of setting mating pool.

**Algorithm 2: Set Mating Pool**

Input: The initialized population P which a set of solution and the current problem Prob 

Output: The mating pool M which a set of solution

```
\1.   Define mating pool M with size of POP_SIZE 

\2.     for i = 0 to POP_SIZE - 1 do

\3.       //calculate the whole fittness of the population

\4.          for j = 0 to POP_SIZE – 1 do

\5.           total_ob += P[i].objective

\6.       end for

\7.       //create the corona

\8.       corona[POP_SIZE]

\9.       for j = 1 to POP_SIZE do

\10.         total_pbb += P[j].objective/total_ob;

\11.        corona[j] = total_pbb;

\12.      end for

\13.      //select the parent

\14.      generate random probability rand_pbb from 0 to 1

\15.      for j = 1 to POP_SIZE do

\16.         if (corona[j] >= rand_pbb)

\17.           index = j-1

\18.           break

\19.         end if

\20.      end for

\21.      M[i] = P[index]

\22.    end for

\23.  return M
```

**3.**   **Mate**

Firstly, I randomly select two parent generations and randomly generate a float number rand_cro from 0 to 1. If rand_cro is less than 0.5, genes are provided by parent1, and conversely by parent2. At the same time there may be a lower probability of mutation during mating. In my program two parents produce one child, the child replaces parent which has the lowest objective in the population. Algorithm 3 is the pseudo code of mating.

**Algorithm 3: Mate**

Input: The mating pool M and the population P which are set of solution

Output: The population which is replaced by child P_new

```
\1.   Define new population P_new with size of POP_SIZE

\2.   genenum = P[0].prob->n

\3.   generate random index rand_par1 from 0 to POP_SIZE-1

\4.   generate random index rand_par2 from 0 to POP_SIZE-1

\5.   parent_1 = M [rand_par1]

\6.   parent_2 = M [rand_par2]

\7.   generate solution child

\8.   //using one-child uniform crossover and mutation

\9.   sort P by objective from low to high

\10.  for k = 0 to popsize

\11.    for i = 0 to genenum

\12.      generate random probability rand_cross from 0 to 1

\13.      generate random probability rand_mut from 0 to 1

\14.       if(rand_cro < 0.5)

\15.         child.x[i] = parent_1.x[i]

\16.       else

\17.         child.x[i] = parent_2.x[i]

\18.      end if

\19.      if(rand_mut <= MUTATION_RATE) // mutation

\20.         if(child.x[i] == 0){

\21.              child.x[i] = 1

\22.         end if

\23.       else

\24.            child.x[i] = 0

\25.       end if

\26.     end for

\27.    evaluate_solution(child)

\28.    if (child.feasibility != 1)

\29.       child = Repair feasibility (child, prob)

\30.    end if

\31.    P[k] = child

\32.  end for

\33.  P_new = P

\34.  return P_new


```

**4.**   **Repair feasibility**

After mating, there might be some solutions exceed the maximum size of knapsack. I use a repair feasibility algorithm to solve this problem. Algorithm4 is the pseudo code of repair feasibility.

**Algorithm 4:** **Repair feasibility**

Input: The child C which needs to be repaired and the current problem Prob

Output: The fixed child C_fix

```
\1.   genenum = Prob->n

\2.   dimnum = Prob->dim

\3.   for i=0 to genenum

\4.     avg_size=0

\5.     item_i = prob->items[i]

\6.     for d=0 to dimnum

\7.       avg_size += item_i->size[d]/prob->capacities[d]

\8.     end for

\9.     item_i->ratio = item_i->p/avg_size

\10.  end for

\11.  sort (items) by ratio from low to high

\12.  for i = 0 to genenum

\13.     if(C.x[prob->items[i].indx] == 1)

\14.        C.x[prob->items[i].indx] = 0

\15.        evaluate_solution(C)

\16.      end if

\17.      if (C.feasibility == 1)

\18.        break

\19.      end if

\20.  end for

\21.  sort (items) by index from low to high

\22.  C_fix = C

\23.  return C_fix;
```



**5.**   **Local search**

I use the Variable Neighborhood Search combine with the best descent algorithm. Comparing with Variable Neighborhood Search combine with the first descent algorithm, it would provide lower gap of my solution.

**Algorithm 5: VNS with best descent algorithm**

Input: initial solution S

Output: best_S

```
\1.   best_S = S

\2.   while (i <=k) 

\3.     S’=best_descent(S); //local search 

\4.     if(S’.objective < S.objective)

\5.       S=S’

\6.       i=1

\7.       continue

\8.     else i++

\9.     endif

\10.  end while

\11.  if S. objective < best_S. objective

\12.    best_S = Send

\13.  end if

\14.  return best_S
```

**6.**   **Memetic** 

In my program, local search is conducted every 200 generations of genetic algorithm, and the Generation gap of each Generation is adjusted by controlling the replace rate. Algorithm6 is the pseudo code of Memetic algorithm.

**Algorithm 6: Memetic algorithm**

Input: The current problem Prob

Output: The best solution S

```
\1.   pop = Algorithm1(Prob)

\2.   time_start = get time

\3.   for i = 0 to 500

\4.     time_fin= get time

\5.     time_spent = time_fin-time_start

\6.     if(time_spent >= MAX_TIME)

\7.       goto end

\8.     end if

\9.     for k = 0 to 200

\10.      mating_pool = Algorithm2(pop, my_pop,prob)

\11.      child = Algorithm3(mating_pool)

\12.      sort (pop) by objective from low to high

 

\13.      //Generation gap

\14.      for j=0 to POP_SIZE*REPLACE_RATE

\15.        my_pop[j] = child 

\16.      end for

\17.    end for

\18.    for l = 0 to POP_SIZE

\19.      Algorithm5 (pop[l])

\20.    end for

\21.  end for     

\22.  end:

\23.  sort (pop) by objective from high to low

\24.  return pop[0]
```

**2.**  **Parameter tuning process**

All test results are obtained under the window operating system, the CPU model is Intel I7-7700HQ. I selected mknapcb5.txt as the test set. When I change a variable, I need to adjust other variables at the same time to reduce its gap. 

**1.**   **POP_SIZE**

With the increase of POP_SIZE, the search breadth of the genetic algorithm will increase rapidly, but at the same time the local search time will also become longer, resulting in a reduction in the number of populations. We need to consider both the breadth of the search and the number of populations to find the best parameters.

| POP_SIZE | GAP(%) |
| -------- | ------ |
| 75       | 3.21   |
| 50       | 2.38   |
| 20       | 1.95   |
| 10       | 1.57   |

 

**2.**   **REPLACE RATE and MUTATION RATE**



A low mutation rate will cause insufficient disturbance, and an excessively high mutation rate will make it approximate to a random search in the late iteration. Low REPLACE_RATE rate will cause the population to evolve slowly, and high REPLACE_RATE will cause excellent individuals in the population to be eliminated.

| REPLACE  RATE | GAP(%) |
| ------------- | ------ |
| 0.8           | 1.46   |
| 0.6           | 2.03   |
| 0.4           | 2.11   |
| 0.2           | 2.13   |

| MUTATION  RATE | GAP(%) |
| -------------- | ------ |
| 0.01           | 2.77   |
| 0.007          | 2.65   |
| 0.005          | 2.44   |
| 0.003          | 1.46   |

 

 

 

 

 

 

 

**3.**   **Different kind of local search**

Here I experiment with two different local searches, VNS with best descent and VNS with first descent. Each of these two local search methods has advantages and disadvantages. VNS with first descent has a faster calculation speed and can quickly generate new populations. Although VNS with best descent takes a long time, it can get more accurate results. We need to evaluate the performance at different times when using two local searches.

|                  | GAP  of VNS with first descent(%) | GAP  of VNS with best descent(%) |
| ---------------- | --------------------------------- | -------------------------------- |
| Max  Time = 100s | 2.37                              | 2.98                             |
| Max  Time = 200s | 2.26                              | 2.43                             |
| Max  Time =300s  | 2.10                              | 1.87                             |
| Max  Time =500s  | 1.97                              | 1.64                             |

 

After the above experiments, when MAXTIME = 300, it is better to choose the VNS with best descent algorithm. And REPLACE RATE should be equal to 0.8 MUTATION RATE equal to 0.003, POP_SIZE equal to 10 can get the best performance of the algorithm.

**3.**  **Result**

Figure1 shows the detailed results of mknapcb1.txt, mknapcb5.txt, mknapcb8.txt.

**![img](file:////Users/dongzhanxie/Library/Group%20Containers/UBF8T346G9.Office/TemporaryItems/msohtmlclip/clip_image001.png)**

figure 1 The best result

The following table shows the best GAP value, the worst GAP value, and the average GAP value of mknapcb2.txt, mknapcb3.txt, mknapcb4.txt, mknapcb6.txt, mknapcb7.txt. 

|                            | mknapcb2 | Mknapcb3 | Mknapcb4 | Mknapcb6 | Mknapcb7 |
| -------------------------- | -------- | -------- | -------- | -------- | -------- |
| the  best GAP value (%)    | 0.01     | 0.04     | 0.23     | 0.56     | 0.76     |
| the  worst GAP value (%)   | 0.96     | 0.86     | 2.02     | 3.02     | 3.92     |
| the  average GAP value (%) | 0.33     | 0.44     | 0.88     | 1.33     | 1.44     |

 

**4.**  **Reflection**

Advantages:

\1.   The memetic algorithm can get many different solutions instead of a single solution.

\2.   Memetic algorithm uses a global method, which greatly expands the search range.

\3.   Strong fault tolerance, even if the difference between the initial population and the optimal solution is very large, after the evolutionary operation and local search, the poorly performing individuals can be eliminated.

\4.   Memetic has strong parallelism, randomness and uncertainty. It can also be compared with other answers in the population to avoid local optimal solutions.

Disadvantages:

\1.   The process of adjusting parameters is very difficult, and we can only rely on step-by-step testing to obtain the best parameters.

\2.   The test results vary greatly on different computer.

\3.   Memetic requires a lot of memory and CPU performance.