#ifndef __PY_ELASTIC_FORCE_FEM_H
#define __PY_ELASTIC_FORCE_FEM_H

#include <elasticForceFEM.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class pyElasticForceFEM : public ElasticForceFEM
{
    public:
    pyElasticForceFEM(py::array_t<float> x, py::array_t<int> tet)
    {
        py::buffer_info buf1 = x.request(), buf2 = tet.request();

        if(buf1.ndim !=1 || buf2.ndim !=1)
        {
            throw std::runtime_error("Number of dimensions must be one");
        }
        
        int xNum = buf1.size / 3;
        int tetNum = buf2.size / 4;

        float *ptr1 = (float*) buf1.ptr;
        int *ptr2 = (int*) buf2.ptr;

        ElasticForceFEM::init(ptr1, xNum, ptr2, tetNum);
    }

    public:
    // setter
    void setDensity(float den)
    {
        this->density = den;
    }

    void setMu01(float _mu01)
    {
        this->mu01 = _mu01;
    }

    void setMu10(float _mu10)
    {
        this->mu10 = _mu10;
    }

    void setV1(float _v1)
    {
        this->v1 = _v1;
    }

    float getDensity()
    {
        return this->density;
    }

    float getMu01()
    {
        return this->mu01;
    }

    float getMu10()
    {
        return this->mu10;
    }

    float gettV1()
    {
        return this->v1;
    }

    float* ComputeFroces(float* u, bool addGravity = false )
    {
        
    }
};

#endif