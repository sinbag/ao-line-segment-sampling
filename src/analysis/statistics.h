//
//  statistics.h
//  PBRT-V3
//
//  Created by Gurprit Singh on 10/6/16.
//
//

#ifndef statistics_h
#define statistics_h

inline double compute_mean(double &mean, const double &integral, int trial)
{
    if(trial == 0){
        std::cerr <<"Progressive mean trial index must start from 1, and not 0 !!!" << std::endl;
        exit(-2);
    }
    else{
        double delta = integral - mean;
        mean += delta / trial;
        //mean = ((trial-1)/double(trial)) * mean + (1/double(trial)) * integral;
    }
    
    return mean;
}

inline double compute_variance(double &variance, const double &mean, const double &integral, int trial)
{
    if(trial == 0){
        std::cerr <<"Progressive variance trial index must start from 1, and not 0 !!!" << std::endl;
        exit(-2);
    }
    else if(trial < 2){
        variance = 0;
    }
    else{
        double deviation = integral - mean;
        variance = ((trial-1)/double(trial)) * variance + (1/double(trial-1)) * deviation * deviation;
    }
    return variance;
}

///
/// The following approach is much less prone to catastrophic cancellation
/// source: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
///
inline void compute_mean_variance(double &mean, double &variance, const double &integral, int trial){
    
    if(trial == 0){
        std::cerr <<"Progressive mean/variance trial index must start from 1, and not 0 !!!" << std::endl;
        exit(-2);
    }
    double delta = integral - mean;
    mean += (delta / trial);
    variance += (delta * (integral - mean));
}


#endif /* statistics_h */
