#pragma once
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>



int config_mapping(int argc, char** argv, boost::program_options::variables_map& dictionary);
int cmdline_settings(int argc, char** argv, double* radius, double* thetamax,
                    double* a, double* b, int* N);


class Configuration {
    public:
        // Public Constructors/Destructors/Functions
        Configuration(int argc, char** argv, std::string file):file(file) {
            std::cout << "construct\n";
            config_mapping(argc,argv,vm);
        };

        int cmdline_settings(int argc, char** argv, double* radius, double* thetamax,
                    double* a, double* b, int* N);
        int config_mapping(int argc, char** argv, boost::program_options::variables_map& dictionary);

        // Public variables
        boost::program_options::variables_map vm;
        std::string file;
        double r = vm["init.conds.radius"].as<double>();
        double dr = vm["init.conds.dr"].as<double>();
        double tmax = vm["init.conds.thetamax"].as<double>();
        double xcenter = vm["init.conds.xcenter"].as<double>();
        double ycenter = vm["init.conds.ycenter"].as<double>();
        int numinit = vm["init.conds.numinit"].as<int>();
        std::string sysfile = vm["system.file"].as<std::string>();
        std::string sysmethod = vm["system.method"].as<std::string>();
        double maxtime = vm["system.timemax"].as<double>();
        double mintime = vm["system.timemin"].as<double>();
        double dt = vm["system.dt"].as<double>();
        float bkg_r = vm["gui.bkg_r"].as<float>();
        float bkg_g = vm["gui.bkg_g"].as<float>();
        float bkg_b = vm["gui.bkg_b"].as<float>();
        float bkg_alpha = vm["gui.bkg_alpha"].as<float>();
        float gui_fovy = vm["gui.fovy"].as<float>();
        float gui_aspect = vm["gui.aspect"].as<float>();
        float gui_znear = vm["gui.znear"].as<float>();
        float gui_zfar = vm["gui.zfar"].as<float>();
        float gui_camera_dist = vm["gui.camera_dist"].as<float>();




};
