/*
 *
 * Copyright (C) 2018-  Yewang Chen<ywchen@hqu.edu.cn; 693849393@qq.com; peterhu@hqu.edu.cn>
 * License: GPL v1
 * This software  may be modified and distributed under the terms
 * of license, and is  prohibitted in any commercial use.
 *
 */

#include "cover_tree.h"
#include <windows.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <vector>
#include "cyw_timer.h"
#define MAXFLOAT 1000000

// Compute the k nearest neighbors
typedef struct Point_K{
    float p,q;
    float d;//d is distance between p and q
};

struct node_p{
    int index;
    float value;
};

using namespace std;

int get_dim(char* s, char* delims){
    char *val_str = NULL;
    val_str = strtok( s, delims );
    int dim=0;
    while( val_str != NULL ) {
        dim++;
        val_str = strtok( NULL, delims );
    }
    return dim;
}

float* get_data(char* s, int dim,char* delims){
    float* temp= (float*) malloc (dim*sizeof(float));
    char *val_str = NULL;
    val_str = strtok( s, delims );
    int counter=0;
    while( val_str != NULL ) {
        //printf(val_str);
        temp[counter]=atof(val_str);
        counter++;
        val_str = strtok( NULL, delims );
    }
    return temp;
}

void read_data_dim_size(char* filename, int* data_dim, int* data_size, char* delims){
    int n_size=0;
    int dim=0;

    char s[10000];
    freopen(filename,"r",stdin);
    while(gets(s))
    {
        if (dim==0)
           dim=get_dim(s,delims);
        n_size ++;
    }
    cout << "dim:" << dim << " n_size:" << n_size << "\n";
    *data_dim=dim;
    *data_size=n_size;
    fclose(stdin);
}

float* read_data(char* filename, char* delims){
    int m = 0;
    int dim, n_size;
    read_data_dim_size(filename,&dim, &n_size, delims);

    float* data= (float*) malloc (n_size*dim*sizeof(float));
    freopen(filename,"r",stdin);
    int counter=0;
    char s[10000];
    while(gets(s))
    {
        float* tmp_data= get_data( s, dim,delims);
        memcpy(data+counter*dim,tmp_data,dim*sizeof(float));
        /*
            for (int i=0;i<dim;i++){
               *(data+counter*dim+i)= tmp_data[i];
               //printf("\n%f, ",*(data+counter*dim+i));
            }
        */
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}

float* read_data(char* filename, char* delims, int* dim, int* data_size){
    int m = 0;
    read_data_dim_size(filename,dim, data_size, delims);

    float* data= (float*) malloc ((*data_size)*(*dim)*sizeof(float));
    freopen(filename,"r",stdin);
    int counter=0;
    char s[10000];
    while(gets(s))
    {
        float* tmp_data= get_data( s,*dim,delims);
        memcpy(data+counter*(*dim),tmp_data,(*dim)*sizeof(float));
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}

float* read_data_add_index(char* filename, char* delims, int* dim, int* data_size){
    int m = 0;
    read_data_dim_size(filename,dim, data_size, delims);

    float* data= (float*) malloc ((*data_size)*(*dim +1)*sizeof(float));
    freopen(filename,"r",stdin);
    int counter=0;
    char s[10000];
    while(gets(s))
    {
        float* tmp_data= get_data( s,*dim,delims);
        memcpy(data+counter*(*dim+1),tmp_data,(*dim+1)*sizeof(float));
        *(data+counter*(*dim+1)+*dim) = counter;
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}

vector<int> getRandom(int total){
    srand((int)time(NULL));
    std::vector<int> input = *new std::vector<int>();
    for (int i = 0; i < total; i++) {
        input.push_back(i);
    }
    vector<int> output = *new vector<int>();

    int end = total;
    for (int i = 0; i < total; i++) {
        vector<int>::iterator iter = input.begin();
        int num = rand()%end;
        iter = iter+num;
        output.push_back(*iter);
        input.erase(iter);
        end--;
    }

    return output;
}

float* Generate_data_with_index(float* raw_data,int data_size,int dim,int new_size){
    float* data = (float*)malloc((new_size)*(dim+1)*sizeof(float));
    int counter = 0;
    vector<int> temp(data_size);
    temp = getRandom(data_size);
    for(int i = 0;i < new_size;i++){
        int index = temp[i];
        memcpy(data+counter*(dim+1),raw_data+index*dim,(dim+1)*sizeof(float));
        *(data+counter*(dim+1)+dim) = counter;
        counter++;
    }
    return data;
}

float* generate_random_data(float* source_data,int data_size, int dim, int radom_data_num){
    float* result= new float[radom_data_num*dim];
    for (int i=0;i<radom_data_num;i++){
        int idx= rand();
        int tmp=idx%data_size;
        for (int j=0;j<dim;j++){
            //float noise= rand()/100;
            result[i*dim+j]= source_data[tmp*dim+j] ;
            //result[i*dim+j]=noise;
        }
    }
    return result;
}

float* generate_query_data(float* source_data,int batch,int batch_num,int dim){
    float* query_data = new float[batch_num*dim];
    for(int i = batch*batch_num;i < (batch+1)*batch_num;i++){
        for(int j = 0;j < dim;j++){
            query_data[(i-batch*batch_num)*dim+j] = source_data[i*dim+j];
        }
    }
    cout << "Generate query data:" << batch*batch_num << "--" << (batch+1)*batch_num << endl;
    return query_data;
}

int write_file_log(const char* buffer,  FILE* fp){
    const int b_size=strlen(buffer);
    if (fp==0) {
        std::cout<<"can't open file\n";
        return 0;
    }
    fseek(fp, 0, SEEK_END);

    fwrite(buffer, b_size, 1, fp);
    return 0;
}

void my_strcat(char* buffer, int val){
    char str_tmp[200];
    stringstream s;
    memset(str_tmp, 0, sizeof(str_tmp));
    //itoa(val,str_tmp,10);
    s << val;
    s >> str_tmp;
    strcat(buffer,str_tmp);
}

bool cmp(struct node_p a, struct node_p b){
    if(a.value < b.value){
        return true;
    }
    return false;
}

void quickSort(struct Point_K s[], int l, int r){
    if (l< r)
    {
        int i = l, j = r;
        struct Point_K x = s[l];
        while (i < j)
        {
            while(i < j && s[j].d >= x.d) // 从右向左找第一个小于x的数
                j--;
            if(i < j)
                s[i++] = s[j];
            while(i < j && s[i].d < x.d) // 从左向右找第一个大于等于x的数
                i++;
            if(i < j)
                s[j--] = s[i];
        }
        s[i] = x;
        quickSort(s, l, i - 1); // 递归调用
        quickSort(s, i + 1, r);
    }
}

void Point_Sort(float* data, int K){
    struct Point_K p[K];
    for(int i = 0;i < K;i++){
        p[i].p = data[i*3];
        p[i].q = data[i*3+1];
        p[i].d = data[i*3+2];
    }
    quickSort(p,0,K-1);
    for(int i = 0;i < K;i++){
        data[i*3] = p[i].p;
        data[i*3+1] = p[i].q;
        data[i*3+2] = p[i].d;
    }
}

float distance(float* p1, float* p2, int data_dim){
    float sum = 0.;
    float* end = p1 + data_dim;
    for (; p1 != end; p1++,p2++){
        float d1 = *p1 - *p2;
        d1 *= d1;
        sum = sum + d1;
    }
    return sqrt(sum);
}

void ComputeDistance(v_array<v_array<float> >* result,int dim,int batch,int batch_num,int K,float* distance_re,float* raw_data){
    int k = 0,p_d = dim -1,p1_index,p2_index;//p_d is the dim of point except index
    float d;//d is the distance of p1 between p2
    float *p1, *p2;
    for(int i = 0;i < batch;i++){
        for(int j = 0;j < result[i].index;j++){
            p1_index = (int)result[i][j][0];
            p1 = raw_data + dim*p1_index;
            for(int k = 1;k < result[i][j].index&&k <= K;k++){
                p2_index = (int)result[i][j][k];
                p2 = raw_data + dim*p2_index;
                d = distance(p1,p2,p_d);
                distance_re[3*K*p1_index+3*(k-1)] = p1_index;
                distance_re[3*K*p1_index+3*(k-1)+1] = p2_index;
                distance_re[3*K*p1_index+3*(k-1)+2] = d;
            }
        }
    }
}

void PretreatDistance(float* dis_matrix,int data_size,int K){
    float* temp = new float[3*K];
    for(int i = 0;i < data_size;i++){
        memcpy(temp,dis_matrix+3*K*i,(3*K)*sizeof(float));
        Point_Sort(temp,K);
        memcpy(dis_matrix+3*K*i,temp,(3*K)*sizeof(float));
    }
}

float* ComputeDensity(float* dis_matrix,int data_size,int K){
    float* density = new float[data_size];
    for(int i = 0;i < data_size;i++){
        density[i] = 1.0/dis_matrix[3*(i+1)*K-1];
    }
    return density;
}

int* FindCluster(float* density,float* delta,int cl,int data_size){
    float* r = (float*)malloc(data_size*sizeof(float));
    for(int i = 0;i < data_size;i++)
        r[i] = density[i]*delta[i];
    int* cluster = (int*)malloc(cl*sizeof(int));
    for(int i = 0;i < cl;i++){
        cluster[i] = 0;
        for(int j = 1;j < data_size;j++){
            int flag = -1,cl = cluster[i];
            if(r[j] > r[cl]){
                for(int k = 0;k < i;k++){
                    if(j == cluster[k])  flag = 0;
                }
                if(flag == -1)
                    cluster[i] = j;
            }
        }
    }
    return cluster;
}

void Find_tem_core(float* dis_matrix,float* density,int* tem_Core,float* delta,int data_size,int K){
    for(int i = 0;i < data_size;i++){
        int p1 = dis_matrix[i*K*3];
        for(int j = 1;j < K;j++){
            int p2 = dis_matrix[i*K*3+j*3+1];
            if(density[p2] > density[p1] && p2 < data_size){
                tem_Core[i] = p2;
                delta[i] = dis_matrix[i*K*3+j*3+2];
                break;
            }
        }
    }
}

long PreProcess_local_density_peak(float* density,float* raw_data,float* dis_matrix,float* delta,int* tem_Core,
                            int local_peak_threshold,int data_size,int K,int dim,node node_data){
    int p_d = dim -1;//p_d is the dim of point except index
    int local_peak_num = 0;
    long c_d_n = 0;//c_d_n is the number of computing distance
    float *p1, *p2, d;
    //count the number of local density peak
    for(int i = 0;i < data_size;i++){
        if(tem_Core[i] == -1)
            local_peak_num++;
    }
    if(local_peak_num > local_peak_threshold){
        int fre = 1;//fre is the frequency of loop
        while(local_peak_num > local_peak_threshold){
            std::cout << "local_peak_num:" << local_peak_num << endl;
            int* peak_index = (int *)malloc(local_peak_num*sizeof(int));
            int num = 0;
            //array peak_index stores the index of local density peaks
            for(int i = 0;i < data_size;i++){
                if(tem_Core[i] == -1)
                    peak_index[num++] = i;
            }
            //array peak_data stores the data of local density peaks
            float* peak_data = (float *)malloc(local_peak_num*dim*sizeof(float));
            for(int i = 0;i < local_peak_num;i++)
                memcpy(peak_data+i*dim,raw_data+peak_index[i]*dim,dim*sizeof(float));
            //building tree for local density peaks
            v_array<point> peak_queries = parse_points(peak_data,local_peak_num,dim);
            node node_query = batch_create(peak_queries);

            v_array<v_array<float> > result;
            k_nearest_neighbor(node_data,node_query,result,K*(fre*2),dim);

            int p1_index,p2_index,K_dis_threshold,flag;
            for(int i = 0;i < result.index;i++){
                p1_index = (int)result[i][0];
                p1 = raw_data + p1_index*dim;
                K_dis_threshold = dis_matrix[(p1_index+1)*3*K-1];
                flag = 0;
                for(int j = 1;j < result[i].index&&j <= (fre+1)*K;j++){
                    p2_index = (int)result[i][j];
                    p2 = raw_data + p2_index*dim;
                    d = distance(p1,p2,p_d);
                    c_d_n++;
                    if(d > K_dis_threshold&&density[p2_index] > density[p1_index]){
                        if(flag == 0){
                            tem_Core[p1_index] = p2_index;
                            delta[p1_index] = d;
                            flag = 1;
                        }else{
                            if(d < delta[p1_index]){
                                tem_Core[p1_index] = p2_index;
                                delta[p1_index] = d;
                            }
                        }
                    }
                }
                if(tem_Core[p1_index] != -1)
                    local_peak_num--;
            }
            fre++;
            free(peak_data);
            free(peak_index);
        }
    }
    return c_d_n;
}

int Find_local_density_peak(float* density,int* tem_Core,float* delta,float* raw_data,int local_peak_threshold,long long* C_D_N,
                            float* dis_matrix,int data_size,int dim,int K,char* buffer_l,node node_data, node_p* node_p_ptr){
    CYW_TIMER c1_timer, c2_timer, c3_timer, c4_timer, c5_timer, c6_timer;
    c1_timer.start_my_timer();
    node_p * p = node_p_ptr;
    for(int i = 0;i < data_size;i++){
        p[i].index = i;
        p[i].value = density[i];
    }
    sort(p,p+data_size,cmp);
    c1_timer.stop_my_timer();
    printf("Runtime of sorting density for source data is %f s\n",c1_timer.get_my_timer());
    strcat(buffer_l,"\n|    Sorting density for source data       |            ");
    c1_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    c6_timer.start_my_timer();
    long c_d_n = PreProcess_local_density_peak(density,raw_data,dis_matrix,delta,tem_Core,local_peak_threshold,data_size,K,dim,node_data);
    c6_timer.stop_my_timer();
    printf("Runtime of PreProcessing local density peak by local D_peak threshold is %f s\n",c6_timer.get_my_timer());
    strcat(buffer_l,"\n|    Preprocessing LDP by LPT              |            ");
    c6_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    int p_d = dim -1;//p_d is the dim of point except index
    int local_density_peak_num = 0;
    long long compute_distance_num = c_d_n;
    float *p1, *p2, den;
    int* flag = (int *)malloc(data_size*sizeof(int));
    c2_timer.start_my_timer();
    for(int i = 0;i < data_size;i++){
        if(tem_Core[i] == -1){
            memset(flag,0,data_size*sizeof(int));
            local_density_peak_num++;
            int den_order = 0;
            c5_timer.start_my_timer();
            //find the order of point i by density
            for(int j = 0;j < data_size;j++){
                if(i == p[j].index){
                    den_order = j;
                }
            }
            c5_timer.stop_my_timer();
            c3_timer.start_my_timer();
            //if point den_order of density is the largest than other points
            if(den_order == data_size - 1){
                float max_dis = 0;
                for(int j = 0;j < data_size;j++){
                    p1 = raw_data+i*dim;
                    p2 = raw_data+j*dim;
                    float d = distance(p1,p2,p_d);
                    if(max_dis < d)
                        max_dis = d;
                }
                delta[i] = max_dis;
                c3_timer.stop_my_timer();
                printf("Runtime of finding delta for root is %f s\n",c3_timer.get_my_timer());
                continue;
            }
            c4_timer.start_my_timer();
            p1 = raw_data + (p[den_order].index)*dim;
            float local_nearest_distance = 1000000;
            int local_nearest_index = -1;
            for(int k = den_order + 1;k < data_size;k++){
                int local_index = p[k].index;
                //if local point is not be filtered
                if(flag[local_index] == 0&&density[i] <= density[local_index]){
                    p2 = raw_data + local_index*dim;
                    float d = distance(p1,p2,p_d);
                    compute_distance_num++;
                    //if distance between local point and local density peak is smaller than local nearest distance
                    if(d < local_nearest_distance){
                        local_nearest_distance = d;
                        local_nearest_index = local_index;
                        continue;
                    }
                    //filter some points by triangle inequality
                    for(int j = K;j > 0;j--){
                        float d_p_j = dis_matrix[local_index*3*K+3*j-1];
                        if(d > local_nearest_distance + d_p_j){
                            for(int s = 1;s <= j;s++){
                                int index = dis_matrix[local_index*3*K+3*s-2];
                                flag[index] = 1;
                            }
                            break;
                        }
                    }
                }
            }
            tem_Core[i] = local_nearest_index;
            delta[i] = local_nearest_distance;
            c4_timer.stop_my_timer();
        }
    }
    printf("Runtime of finding order end,runtime is %f s\n",c5_timer.get_my_timer());
    //printf("\nRuntime of finding parent node for local density peak is %f s\n",c4_timer.get_my_timer());
    //std::cout << "compute distance num:" << compute_distance_num << endl;
    //strcat(buffer_l,"\ncompute distance num: ");
    //my_strcat(buffer_l,compute_distance_num);
    *C_D_N = compute_distance_num;
    //std::cout << "local density peak num:" << local_density_peak_num << endl;
    //strcat(buffer_l,"\nlocal density peak num: ");
    //my_strcat(buffer_l,local_density_peak_num);
    //LDP_N = local_density_peak_num;
    c2_timer.stop_my_timer();
    printf("Runtime of finding parent for local density peak is %f s\n",c2_timer.get_my_timer());
    strcat(buffer_l,"\n|    Finding parent nodes for LDP          |            ");
    c2_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    free(flag);
    free(p);
    return local_density_peak_num;
}

void Label_cluster(int* Point_cl,int* tem_Core,int* cluster,int cl,int data_size){
    int c = 0;
    for(int i = 0;i < data_size;i++){
        int flag = -1;//Not find Core in cluster
        int p = i;
        while(flag == -1){
            for(int j = 0;j < cl;j++){
                if(tem_Core[p] == cluster[j]||p == cluster[j]){
                    flag = 0;//Find Core in cluster
                    c = j+1;
                }
            }
            if(flag == 0)
                Point_cl[i] = c;
            else
                p = tem_Core[p];
        }
    }
}

void Save_Result(int* result,int data_size,FILE* file){
    char buffer_r[10];
    for(int i = 0;i < data_size;i++){
        memset(buffer_r,0,sizeof(buffer_r));
        my_strcat(buffer_r,result[i]);
        strcat(buffer_r,"\n");
        write_file_log(buffer_r,file);
    }
    fclose(file);
}

void Fast_Density_Peak(char* data_file_name,int K, float* raw_data,int data_size,int batch_num,int dim,int local_peak_threshold,
                    FILE* log_file, FILE* results,int cl,FILE* cl_results,float* dis_matrix, node_p *node_p_ptr){
    char buffer_l[10000];
    memset(buffer_l,0,sizeof(buffer_l));

    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n|                            FastDPeak Version 1.0                            |");
    strcat(buffer_l,"\n|===============================+=============================================|");
    strcat(buffer_l,"\n|           File name           |               ");
    strcat(buffer_l,data_file_name);
    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n|          Cardinality          |               ");
    my_strcat(buffer_l,data_size);
    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n|           Dimension           |               ");
    my_strcat(buffer_l,dim);
    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n|     Local Peak Threshold      |               ");
    my_strcat(buffer_l,local_peak_threshold);
    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n|               K               |               ");
    my_strcat(buffer_l,K);
    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");

    dim++;
    v_array<point> data_set = parse_points(raw_data,data_size,dim);
    CYW_TIMER build_source_timer, build_query_timer, query_timer, cl_timer, ctree_time;
    CYW_TIMER c1_dis_timer, c2_dis_timer, c3_dis_timer, c4_dis_timer, c5_dis_timer;

    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n| Processes:                                             Running              |");
    strcat(buffer_l,"\n|    Process name                                         Time(s)             |");
    strcat(buffer_l,"\n+==========================================+==================================+");

    printf("building tree for source data....\n");
    build_source_timer.start_my_timer();
    node node_data = batch_create(data_set);
    build_source_timer.stop_my_timer();
    printf("Runtime of building tree for source data is %f s\n",build_source_timer.get_my_timer());
    strcat(buffer_l,"\n|    Building tree for source data         |            ");
    build_source_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    int batch = data_size/batch_num;
    std::cout << "batch:" << batch << endl;

    cl_timer.start_my_timer();
    v_array<v_array<float> >* res = new v_array<v_array<float> >[batch];
    for (int j=0;j<batch;j++){
        v_array<v_array<float> >* tmp_batch= new v_array<v_array<float> >(batch_num);
        res[j]= *tmp_batch;
        for (int i=0;i<batch_num;i++){
            v_array<float> tmp_v(K);
            push(res[j],tmp_v);
        }
        res[j].index=0;
    }

    ctree_time.start_my_timer();
    for(int i = 0;i < batch;i++){
        std::cout<<"generate queries...\n";
        float* query_data = generate_query_data(raw_data,i,batch_num,dim);
        v_array<point> queries = parse_points(query_data,batch_num,dim);

        printf("building tree for query data....\n");
        build_query_timer.start_my_timer();
        node node_query = batch_create(queries);
        build_query_timer.stop_my_timer();
        printf("Runtime of building tree for query data is %f s\n",build_query_timer.get_my_timer());
        if(i == batch-1){
            strcat(buffer_l,"\n|    Building tree for query data          |            ");
            build_query_timer.strcat_to_buffer(buffer_l);
        }
        //strcat(buffer_l,"s");

        query_timer.start_my_timer();
        k_nearest_neighbor_new(node_data,node_query,res[i],K,dim);
        query_timer.stop_my_timer();
        printf("K = %d, runtime of kNN is %f s\n",K,query_timer.get_my_timer());

        /*
        strcat(buffer_l,"\n K= ");
        my_strcat(buffer_l, K);
        strcat(buffer_l,", runtime of kNN:");
        query_timer.strcat_to_buffer(buffer_l);
        strcat(buffer_l,"s");
        */
        queries.free_resource();
    }
    data_set.free_resource();
    ctree_time.stop_my_timer();

    printf("K = %d,runtime of building covertree and KNN is %f s\n",K,ctree_time.get_my_timer());
    strcat(buffer_l,"\n|    Building covertree and KNN            |            ");
    ctree_time.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    c1_dis_timer.start_my_timer();
    ComputeDistance(res,dim,batch,batch_num,K,dis_matrix,raw_data);
    c1_dis_timer.stop_my_timer();
    printf("Runtime of computing distance for query data: %f s\n",c1_dis_timer.get_my_timer());
    strcat(buffer_l,"\n|    Computing distance for query data     |            ");
    c1_dis_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    for (int i = 0;i < res->length;i++){
        res->elements[i].free_resource();
    }
    res->free_resource();

    int d_size = batch * batch_num;
    PretreatDistance(dis_matrix,d_size,K);

    float* density = ComputeDensity(dis_matrix,d_size,K);

    int* tem_core = (int *)malloc(d_size*sizeof(int));
    float* delta = (float *)malloc(d_size*sizeof(float));
    memset(tem_core,-1,d_size*sizeof(int));
    memset(delta,-1,d_size*sizeof(float));

    c2_dis_timer.start_my_timer();
    Find_tem_core(dis_matrix,density,tem_core,delta,d_size,K);
    c2_dis_timer.stop_my_timer();
    printf("Runtime of finding initiative LDP for data %f s\n",c2_dis_timer.get_my_timer());
    strcat(buffer_l,"\n|    Finding initiative LDP for data       |            ");
    c2_dis_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");

    long long C_D_N = 0;
    int LDP_N = 0;

    LDP_N = Find_local_density_peak(density,tem_core,delta,raw_data,local_peak_threshold,&C_D_N,dis_matrix,d_size,dim,K,buffer_l,node_data,node_p_ptr);
    //Find_local_density_peak(density,tem_core,delta,raw_data,local_peak_threshold,dis_matrix,d_size,dim,K,buffer_l,node_data, node_p_ptr);

    c3_dis_timer.start_my_timer();
    int* cluster = FindCluster(density,delta,cl,d_size);
    /*
    for(int i = 0;i < cl;i++)
        std::cout << i+1 << ":" << cluster[i] << endl;
    */
    c3_dis_timer.stop_my_timer();
    printf("Runtime of determining final clusters %f s\n",c3_dis_timer.get_my_timer());
    strcat(buffer_l,"\n|    Determing final clusters              |            ");
    c3_dis_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");
    delete density;
    free(delta);

    c4_dis_timer.start_my_timer();
    int* Point_cl = (int*)malloc(data_size*sizeof(int));
    Label_cluster(Point_cl,tem_core,cluster,cl,data_size);
    c4_dis_timer.stop_my_timer();
    printf("Runtime of labeling cluster for data is %f s\n",c4_dis_timer.get_my_timer());
    strcat(buffer_l,"\n|    Labeling cluster for data             |            ");
    c4_dis_timer.strcat_to_buffer(buffer_l);
    //strcat(buffer_l,"s");
    strcat(buffer_l,"\n+------------------------------------------+----------------------------------+");

    cl_timer.stop_my_timer();
    printf("Runtime of clustering for data is %f s\n",cl_timer.get_my_timer());

    strcat(buffer_l,"\n+-----------------------------------------------------------------------------+");
    strcat(buffer_l,"\n|                          Result of FastDPeak clustering                     |");
    strcat(buffer_l,"\n+===============================+=============================================+");
    strcat(buffer_l,"\n|        Overall runtime        |                   ");
    cl_timer.strcat_to_buffer(buffer_l);
    strcat(buffer_l,"s");
    strcat(buffer_l,"\n+-------------------------------+---------------------------------------------+");
    strcat(buffer_l,"\n|      Compute distance num     |                     ");
    my_strcat(buffer_l,C_D_N);
    strcat(buffer_l,"\n+-------------------------------+---------------------------------------------+");
    strcat(buffer_l,"\n|            LDP num            |                     ");
    my_strcat(buffer_l,LDP_N);
    strcat(buffer_l,"\n+-------------------------------+---------------------------------------------+");
    strcat(buffer_l,"\n|            cluster                                 Points                   |");
    strcat(buffer_l,"\n+===============================+=============================================+");

    int *cl_points = (int *)malloc((cl+1)*sizeof(int));
    for(int i = 0;i < cl+1;i++){
        cl_points[i] = 0;
    }
    int temp;
    for(int i = 0;i < data_size;i++){
        temp = Point_cl[i];
        cl_points[temp] ++;
    }
    for(int i = 1;i <= cl;i++){
        strcat(buffer_l,"\n|               ");
        my_strcat(buffer_l,i);
        strcat(buffer_l,"               |                     ");
        my_strcat(buffer_l,cl_points[i]);
    }
    strcat(buffer_l,"\n+-------------------------------+---------------------------------------------+");

    free(cl_points);
    free(cluster);

    write_file_log(buffer_l,log_file);
    fclose(log_file);

    Save_Result(tem_core,d_size,results);
    Save_Result(Point_cl,d_size,cl_results);

    free(tem_core);
    free(Point_cl);
}

int main(int argc, char *const argv[]){

    //Input data file
    char *data_file_name = "data/synthesis_3.txt";
    //Output log_files and result files
    FILE* log_file = fopen("exp_log.txt","a+");
    FILE* results= fopen("result.txt","a+");
    FILE* cl_results = fopen("cl_results.txt","a+");

    std::cout << "reading data...\n";
    //data_size is the size of raw_data file,new_size is the size of new_data
    int dim, data_size, new_size = 2000;
    //read data from data_file
    float *raw_data = read_data(data_file_name," ",&dim, &data_size);
    //generate new_data with index by new_size
    float *data_with_index = Generate_data_with_index(raw_data,data_size,dim,new_size);
    free(raw_data);
    std::cout << "data read";

    int K = 11,batch_num = 1000,cl = 5,local_peak_threshold = 30;
    float *dis_matrix = (float*)malloc(3*data_size*K*sizeof(float));
    node_p *node_p_ptr = new node_p [data_size];
    //perform Fast Density Peak Clustering
    Fast_Density_Peak(data_file_name,K,data_with_index,new_size,batch_num,dim,local_peak_threshold
                      ,log_file,results,cl,cl_results,dis_matrix,node_p_ptr);

    free(data_with_index);
    free(dis_matrix);

    return 0;
}

