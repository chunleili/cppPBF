#include <vector>
#include <array>
#include "vec.h"
#include "myPrint.h"

using namespace zeno;

inline float t_kernelPoly6(float dist, float h)
{
    float coeff = 315.0 / 64.0 / 3.14159265358979323846;
    float res = 0.0;
    if(dist > 0 && dist < h)
    {
        float x = (h * h - dist * dist) / (h * h * h);
        res = coeff * x * x * x;
    }
    return res;
}


inline float t_computeScorr(const vec3f& distVec)
{
    float coeffDq = 0.3;
    float coeffK = 0.3;
    float h = 1.1;
    float x = t_kernelPoly6(length(distVec), h) / t_kernelPoly6(coeffDq * h, h);
    x = x * x;
    x = x * x;
    return (-coeffK) * x;
}



inline vec3f t_kernelSpikyGradient(const vec3f& r, float h)
{
    float coeff = -45.0 / 3.14159265358979323846;
    vec3f res{0.0, 0.0, 0.0};
    float dist = length(r);
    if (dist > 0 && dist < h)
    {
        float x = (h - dist) / (h * h * h);
        float factor = coeff * x * x;
        res = r * factor / dist;
    }
    return res;
}


void test(
    const std::vector<std::vector<int>> &neighborList,
    const std::vector<vec3f> &pos,
    std::vector<vec3f> &dpos
)
{
    printVectorField(pos,"input_pos.txt",16);
    printVectorField(neighborList, "input_neighborList.txt");

    int numParticles = 100;
    float h = 1.1;
    float mass = 1.0;
    float rho0 = 1.0;

    std::vector<float> lambda(numParticles);
    //compute lambda
    // lambda.clear();
    // lambda.resize(numParticles);
    float lambdaEpsilon = 100.0; // to prevent the singularity
    for (size_t i = 0; i < numParticles; i++)
    {
        vec3f gradI{0.0, 0.0, 0.0};
        float sumSqr = 0.0;
        float densityCons = 0.0;

        for (size_t j = 0; j < neighborList[i].size(); j++)
        {
            int pj = neighborList[i][j];
            if (pj<0) break;
            vec3f distVec = pos[i] - pos[pj];
            vec3f gradJ = t_kernelSpikyGradient(distVec, h);
            gradI += gradJ;
            sumSqr += dot(gradJ, gradJ);
            densityCons += t_kernelPoly6(length(distVec), h);
        }
        densityCons = (mass * densityCons / rho0) - 1.0;

        //compute lambda
        sumSqr += dot(gradI, gradI);
        lambda[i] = (-densityCons) / (sumSqr + lambdaEpsilon);
    }


    //compute dpos
    // dpos.clear();
    // dpos.resize(numParticles);
    for (size_t i = 0; i < numParticles; i++)
    {
        vec3f dposI{0.0, 0.0, 0.0};
        for (size_t j = 0; j < neighborList[i].size(); j++)
        {
            int pj = neighborList[i][j];
            if (pj<0) break;
            vec3f distVec = pos[i] - pos[pj];

            float sCorr = t_computeScorr(distVec);
            dposI += (lambda[i] + lambda[pj] + sCorr) * t_kernelSpikyGradient(distVec, h);
        }
        dposI /= rho0;
        dpos[i] = dposI;
    }

    printVectorField(dpos, "output_dpos.txt",16);
}



// void PBF::testScorr()
// {
//     float dx = 0.00001;
//     float dy = 0.00001;
//     float dz = 0.00001;
//     std::vector<float> arr;
//     for (size_t i = 0; i < 100; i++)
//     {
//         for (size_t j = 0; j < 100; j++)
//         {
//             for (size_t k = 0; k < 100; k++)
//             {
//                 float x = i * dx;
//                 float y = j * dy;
//                 float z = k * dz;
                                
//                 vec3f val = {x,y,z}; 
//                 arr.push_back(computeScorr(val));
//             }
//         }
//     }
//     printScalarField(arr,"sCorr.txt",16);
// }
