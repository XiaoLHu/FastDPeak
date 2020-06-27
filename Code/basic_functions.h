/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */

#ifndef BASIC_FUNTIONS_INCLUDED
#define BASIC_FUNTIONS_INCLUDED
#include<stdio.h>
#include <stdlib.h>
#include "conio.h"
#include "Array.hh"


template<typename T>
inline FLOAT_TYPE scalar_product(const Vector<T>& x, const Vector<T>& y)
{
  register int i, n = x.size();
  register double sum;

  sum = 0.0;
  for (i = 0; i < n; i++)
    sum += x[i] * y[i];
  return sum;
}

// Utility functions for printing vectors and matrices
template<typename T>
void print_matrix(char* name, const Matrix<T>& A, int n = -1, int m = -1);

template<typename T>
void print_matrix(char* name, const Matrix<T>& A, int n, int m)
{
      std::ostringstream s;
      std::string t;
      if (n == -1)
          n = A.nrows();
      if (m == -1)
          m = A.ncols();

      s << name << ": " << std::endl;
      for (int i = 0; i < n; i++)
      {
         s << " ";
         for (int j = 0; j < m; j++)
           s << A[i][j] << ", ";
         s << std::endl;
      }
      t = s.str();
      t = t.substr(0, t.size() - 3); // To remove the trailing space, comma and newline

      std::cout << t << std::endl;
}

template<typename T>
inline FLOAT_TYPE norm_v(Vector<T>& v){
    register T result=0;
    register int i, n=v.size();
    for ( i=0;i<n;i++)
        result += pow( v[i],2);
    return sqrt(result);
}

template<typename T>
inline FLOAT_TYPE norm_v_squre(Vector<T>& v){
    register T result=0;
    register int i, n=v.size();
    for ( i=0;i<n;i++)
        result += pow( v[i],2);
    return result;
}

template<typename T>
inline FLOAT_TYPE pdist2(Vector<T>& a, Vector<T>&b)
{
    Vector<T> c=a-b;
    double result= norm_v<T>(c);
    return result;
}

template<typename T>
Vector<T> pdist2(Matrix<T>& mat, Vector<T>&b){
    register int i, j, rows=mat.nrows();
    register int cols= mat.ncols();

    if (cols!=b.size())
        throw std::runtime_error("Matrix rows is not consistent with b");
    Vector<T> result(mat.nrows());
    for ( i=0;i<rows;i++){
        result[i]=0;
        for (j=0;j<cols;j++)
        {
            result[i]+=pow(mat[i][j]-b[j],2);
        }
        result[i]=sqrt(result[i]);
    }
    return result;
}

template<typename T>
Vector<T> pdist2(Matrix<T>& mat, T* b, int size){
    register int i, cols= mat.ncols(), rows=mat.nrows();
    if (cols!=size)
        throw std::runtime_error("Matrix cols is not consistent with b");
    Vector<T> result(mat.nrows());
    for ( i=0;i<rows;i++){
        result[i]=0;
        for (int j=0;j<size;j++)
        {
            result[i] += pow(mat[i][j]-b[j],2);
        }
        result[i]= sqrt(result[i]);
    }
    return result;
}



template<typename T>
Vector<T> pdist2(Vector<T>&b, Matrix<T>& mat){
    return pdist2 (mat,b);
}

template<typename T>
Vector<T> pdist2(Vector<T>&b, Matrix<T>& mat, Vector<int>& indexes){
    register int i, cols= mat.ncols(),n=indexes.size();
    if (cols!=b.size())
        throw std::runtime_error("Matrix cols is not consistent with b");
    Vector<T> result(n);
    for ( i=0;i<n;i++){
        result[i]=0;
        for (int j=0;j<cols;j++)
        {
            result[i] += pow(mat[indexes[i]][j]-b[j],2);
        }
        result[i]= sqrt(result[i]);
    }
    return result;
}


template<typename T>
inline FLOAT_TYPE pdist2(T* a, T* b ,int size)
{
    register int i;
    double result;
    for ( i=1;i<size;i++){
        result+= pow((a[i]-b[i]),2);
    }
    return sqrt(result);
}


template<typename T>
inline FLOAT_TYPE pdist2_squre(Vector<T>& a, Vector<T>&b)
{
    Vector<T> c=a-b;
    double result= norm_v_squre<T>(c);
    return result;
}

template<typename T>
void pdist2_squre(Matrix<T>& mat, Vector<T>&b, Vector<T>& result, int from_loc){
    register int i, j, rows=mat.nrows();
    register int cols= mat.ncols();

    if (cols!=b.size())
        throw std::runtime_error("Matrix rows is not consistent with b");
    for ( i=0;i<rows;i++){
        result[from_loc+i]=0;
        for (j=0;j<cols;j++)
        {
            result[from_loc+i]+=pow(mat[i][j]-b[j],2);
        }
    }
}

template<typename T>
Vector<T> pdist2_squre(Matrix<T>& mat, Vector<T>&b){
    register int i, j, rows=mat.nrows();
    register int cols= mat.ncols();

    if (cols!=b.size())
        throw std::runtime_error("Matrix rows is not consistent with b");
    Vector<T> result(mat.nrows());
    for ( i=0;i<rows;i++){
        result[i]=0;
        for (j=0;j<cols;j++)
        {
            result[i]+=pow(mat[i][j]-b[j],2);
        }
    }
    return result;
}

template<typename T>
Vector<T> pdist2_squre(Matrix<T>& mat, T* b, int size){
    register int i, cols= mat.ncols(), rows=mat.nrows();
    if (cols!=size)
        throw std::runtime_error("Matrix cols is not consistent with b");
    Vector<T> result(mat.nrows());
    for ( i=0;i<rows;i++){
        result[i]=0;
        for (int j=0;j<size;j++)
        {
            result[i] += pow(mat[i][j]-b[j],2);
        }
    }
    return result;
}



template<typename T>
Vector<T> pdist2_squre(Vector<T>&b, Matrix<T>& mat){
    return pdist2_squre (mat,b);
}

template<typename T>
Vector<T> pdist2_squre(Vector<T>&b, Matrix<T>& mat, Vector<int>& indexes){
    register int i, cols= mat.ncols(),n=indexes.size();
    if (cols!=b.size())
        throw std::runtime_error("Matrix cols is not consistent with b");
    Vector<T> result(n);
    for ( i=0;i<n;i++){
        result[i]=0;
        for (int j=0;j<cols;j++)
        {
            result[i] += pow(mat[indexes[i]][j]-b[j],2);
        }
    }
    return result;
}


template<typename T>
inline FLOAT_TYPE pdist2_squre(T* a, T* b ,int size)
{
    register int i;
    double result;
    for ( i=1;i<size;i++){
        result+= pow((a[i]-a[b]),2);
    }
}

template<typename T>
void print_vector(char* name, const Vector<T>& v, int n=-1);


template<typename T>
void print_vector(char* name, const Vector<T>& v, int n)
{
  std::ostringstream s;
  std::string t;
  if (n == -1)
    n = v.size();

  s << name << ": " << std::endl << " ";
  for (int i = 0; i < n; i++)
  {
    s << v[i] << ", ";
  }
  t = s.str();
  t = t.substr(0, t.size() - 2); // To remove the trailing space and comma

  std::cout << t << std::endl;
}

template<typename T>
Vector<T> createAnIndexedVector(int start_loc, int end_loc){
    int n=end_loc-start_loc;
    Vector <T> result (n);
    for (int i=0;i<n;i++)
        result[i] =i+start_loc;
}

template<typename T>
Vector<T> createAnIndexedVector(int num){
    Vector <T> result (num);
    for (int i=0;i<num;i++)
        result[i] =i;
    return result;
}


template<typename T>
T g_getBeta(Vector<T>& alpha, Vector<T> &point, Matrix<T>& data){
    T beta= -scalar_product(alpha, point);
    //shift the hyperplane
    beta=g_getBeta_1(alpha,beta,data);
    return beta;
}

template<typename T>
T g_getBeta_1(Vector<T> &alpha, T beta, Matrix<T>& data){
    int rows=data.nrows();
    Vector <T> beta_v (data.nrows());
    Vector<T> dists= alpha.product(data);
    for (int i=0;i<rows;i++)
        dists[i]= dists[i] +beta;

    unsigned min_index= index_min(dists);

    //shift the hyperplane to newPoint
    double newbeta=-scalar_product(alpha, data.extractRow(min_index));
    return newbeta;
}

void my_strcat(char* buffer, int val){
    char str_tmp[200];
    memset(str_tmp, 0, sizeof(str_tmp));
    itoa(val,str_tmp,10);
    strcat(buffer,str_tmp);
}



void my_strcat_longlong(char* buffer, long long val){
    char str_tmp[200];
    sprintf(str_tmp,"%lld",val);

    strcat(buffer,str_tmp);
}

void my_strcat_double(char* buffer, double val){
    char str_tmp[200];
    sprintf(str_tmp,"%f",val);

    strcat(buffer,str_tmp);
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

#endif // BASIC_FUNTIONS_INCLUDED
