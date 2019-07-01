/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */
#ifndef CYW_TIMER_H_INCLUDED
#define CYW_TIMER_H_INCLUDED
void my_strcat_double(char* buffer, double val){
    char str_tmp[200];
    sprintf(str_tmp,"%f",val);

    strcat(buffer,str_tmp);
}

class CYW_TIMER
{
public:
    CYW_TIMER(){
        init_my_timer();
    }
    /**
     * Helper method for computing the current time (w.r.t to an offset).
     *
     *@return System in in microseconds
     */
    long get_system_time_in_microseconds(){
        /* --------RUN IN  LINUX----------
            struct timeval tempo;
            gettimeofday(&tempo, NULL);
            return tempo.tv_sec * 1000000 + tempo.tv_usec;
        */
        LARGE_INTEGER nFreq;
        LARGE_INTEGER t1;
        double dt;
        QueryPerformanceFrequency(&nFreq);
        QueryPerformanceCounter(&t1);
        dt = t1.QuadPart/(double)nFreq.QuadPart;
        return dt* 1000000;
    }

    double get_system_time_in_seconds(){
        /* --------RUN IN  LINUX----------
            struct timeval tempo;
            gettimeofday(&tempo, NULL);
            return tempo.tv_sec * 1000000 + tempo.tv_usec;
        */
        LARGE_INTEGER nFreq;
        LARGE_INTEGER t1;
        double dt;
        QueryPerformanceFrequency(&nFreq);
        QueryPerformanceCounter(&t1);
        dt = t1.QuadPart/(double)nFreq.QuadPart;
        return dt;
    }

    /**
     * Initializes a timer
     *
     * @param *timer Pointer to timer struct instance
     */
    void init_my_timer(){
        start_time = 0.0;
        elapsed_time = 0.0f;
        elapsed_time_total = 0.0f;
    }

    /**
     * Starts a given timer
     *
     * @param *timer Pointer to timer struct instance
     */
    void start_my_timer(){
        start_time = get_system_time_in_seconds();
    }

    /**
     * Stops a given timer
     *
     * @param *timer Pointer to timer struct instance
     */
    void stop_my_timer(){
        elapsed_time = (double)get_system_time_in_seconds() - start_time;
        elapsed_time_total += elapsed_time;
    }

    /**
     * Returns the time measured by a given timer
     *
     * @param *timer Pointer to timer struct instance
     * @return Passed time in seconds
     */
    double get_my_timer(){
        return elapsed_time_total;
    }

    void print(char* out_info){
        std::cout<<out_info;
        print();
    }

    void print(){
        std::cout<<"running time:"<< this->elapsed_time_total*1000000<<"us, "
                <<this->elapsed_time_total*1000<<"ms, " <<this->elapsed_time_total<<"s\n";
    }

    void strcat_to_buffer(char* buffer){
        strcat(buffer,"running time: ");
        my_strcat_double(buffer,elapsed_time_total);
    }

    void strcat_to_buffer(char* out_info, char* buffer){
        strcat(buffer,out_info);
        my_strcat_double(buffer,elapsed_time_total);
    }
private:
    double start_time;
	double elapsed_time;
	double elapsed_time_total;
};


#endif // CYW_TIMER_H_INCLUDED
