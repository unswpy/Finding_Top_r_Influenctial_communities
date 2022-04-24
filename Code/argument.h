//
// Created by sbian on 2020/7/20.
//

#ifndef WEIGHT_COMMUNITY_ARGUMENT_H
#define WEIGHT_COMMUNITY_ARGUMENT_H

class Argument
{
public:
    int _k = 20; // degree constraint k.
    int _r = 20; // top r constraint.
    int _s = 20; // size constraint.
    double _p = 0.0; // local constraint
    double _eps = 0.01; // approximation ratio
    double _tau = 0.02; // threshold constraint tau.
    std::string _func = "sum"; // Aggregated function. sum, avg, min, max
    std::string _graphname = "dblp"; // Graph name. Default is "dblp"
    std::string _mode = "top"; // mode top, enum, local
    std::string _algName = "naive"; // Algorithm. Default is naive. naive, improve, heuristic, approx
    std::string _dir = "spreprocess"; // Directory
    std::string _resultFolder = "result"; // Result folder.
    std::string _outFileName; // File name of the result

    Argument(int argc, char * argv[])
    {
        std::string param, value;
        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') {
                break;
            }
            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');
            if (!param.compare("-k")) {
                _k = stoi(value);
            }
            else if (!param.compare("-r")) {
                _r = stoi(value);
            }
            else if (!param.compare("-s")) {
                _s = stoi(value);
            }
            else if (!param.compare("-p")) {
                _p = stod(value);
            }
            else if (!param.compare("-eps")) {
                _eps = stod(value);
            }
            else if (!param.compare("-tau")) {
                _tau = stod(value);
            }
            else if (!param.compare("-func")) {
                _func = value;
            }
            else if (!param.compare("-gname")) {
                _graphname = value;
            }
            else if (!param.compare("-mode")) {
                _mode = value;
            }
            else if (!param.compare("-alg")) {
                _algName = value;
            }
            else if (!param.compare("-dir")) {
                _dir = value;
            }
            else if (!param.compare("-outpath")) {
                _resultFolder = value;
            }
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;

#endif //WEIGHT_COMMUNITY_ARGUMENT_H
