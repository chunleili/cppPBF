
#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip> 

template<typename T>
void printVector(T contents, int precision=8, int maxTimes=20,  std::string msg="", bool toFile=false, std::string fileName="debugOutput.txt")
{
    static int times = 0;
    times++;
    if (times>maxTimes)
        return;
    
    if(!toFile)
    {
        std::cout<<msg;
        for(auto x:contents)
            std::cout<<std::fixed <<std::setprecision(precision)<<x<<"\t";
        std::cout<<"\n";
        return;
    }
    else if(toFile)
    {
        static std::ofstream fout;
        fout.open(fileName, std::ios::app);
        fout<<msg;
        for(const auto& x:contents)
            fout<<std::fixed <<std::setprecision(precision)<<x<<"\t";
        fout<<"\n";
        fout.close();
        return;
    }
}

template<typename T>
void printVectorField(T content, std::string fileName, size_t precision=8)
{
    std::ofstream f;
    f.open(fileName);
    for(const auto& x:content)
    {
        for(const auto& xx:x)
            f<<std::fixed <<std::setprecision(precision)<<xx<<"\t";
        f<<"\n";
    } 
    f.close();
}

template<typename T>
void printScalarField(T content, std::string fileName, size_t precision=8)
{
    std::ofstream f;
    f.open(fileName);
    for(const auto& x:content)
    {
        f<<std::fixed <<std::setprecision(precision)<<x<<"\t";
        f<<"\n";
    } 
    f.close();
}

#define echo(content) {std::cout<<(#content)<<": "<<content<<std::endl;}