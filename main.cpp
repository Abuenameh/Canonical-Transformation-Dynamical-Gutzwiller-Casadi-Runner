/* 
 * File:   main.cpp
 * Author: Abuenameh
 *
 * Created on 10 January 2015, 17:14
 */

#include <cstdlib>

using namespace std;

#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/date_time.hpp>
#include <boost/process.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace boost;
using namespace boost::filesystem;
using namespace boost::posix_time;
using namespace boost::process;
using namespace boost::process::initializers;
using namespace boost::iostreams;

#include <casadi/casadi.hpp>

#include "mathematica.hpp"

struct arguments {
    int seed;
    double Wi;
    double Wf;
    double mu;
    double D;
    int resi;
};

void combine(ostream& os, string res, int resi, int len) {
    os << res << "[" << resi << "]=Join@@Table[" << res << "[" << resi << ",i],{i,0," << len-1 << "}];" << endl;
}

void threadfunc(int subresi, double taui, double tauf, int ntaus, arguments inargs, path resdir) {
#ifdef AMAZON
        string prog = "/home/ubuntu/canonical_transformat_dynamical_gutzwiller_casadi";
#else
        string prog = "/Users/Abuenameh/NetBeansProjects/Canonical Transformation Dynamical Gutzwiller Casadi/dist/Release/Clang-MacOSX/canonical_transformation_dynamical_gutzwiller_casadi";
#endif
	vector<string> args;
	args.push_back(prog);
        args.push_back(lexical_cast<string>(inargs.seed));
        args.push_back(lexical_cast<string>(inargs.Wi));
        args.push_back(lexical_cast<string>(inargs.Wf));
        args.push_back(lexical_cast<string>(inargs.mu));
        args.push_back(lexical_cast<string>(inargs.D));
        args.push_back(lexical_cast<string>(taui));
        args.push_back(lexical_cast<string>(tauf));
        args.push_back(lexical_cast<string>(ntaus));
        args.push_back(lexical_cast<string>(1));
        args.push_back(lexical_cast<string>(inargs.resi));
        args.push_back(lexical_cast<string>(subresi));
//        args.push_back(lexical_cast<string>());

    ostringstream oss;
    oss << "res." << inargs.resi << "." << subresi << ".out.txt";
    path outfile = resdir / oss.str();
//    filesystem::ofstream os(resfile);
    file_descriptor_sink sink(outfile.string());

//    cout << args << endl;
    child c = execute(set_args(args), inherit_env(), bind_stdout(sink));
        wait_for_exit(c);
}

/*
 * 
 */
int main(int argc, char** argv) {

    ptime begin = microsec_clock::local_time();

    int seed = lexical_cast<int>(argv[1]);
    double Wi = lexical_cast<double>(argv[2]);
    double Wf = lexical_cast<double>(argv[3]);
    double mu = lexical_cast<double>(argv[4]);
    double D = lexical_cast<double>(argv[5]);
    double taui = lexical_cast<double>(argv[6]);
    double tauf = lexical_cast<double>(argv[7]);
    int ntaus = lexical_cast<int>(argv[8]);
    int numthreads = lexical_cast<int>(argv[9]);
    int resi = lexical_cast<int>(argv[10]);
    
    arguments args;
    args.seed = seed;
    args.Wi = Wi;
    args.Wf = Wf;
    args.mu = mu;
    args.D = D;

#ifdef AMAZON
//    path resdir("/home/ubuntu/Results/Canonical Transformation Dynamical Gutzwiller");
    path resdir("/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/Canonical Transformation Dynamical Gutzwiller");
    string loadfunc = "LoadAmazonSubResult";
#else
    path resdir("/Users/Abuenameh/Documents/Simulation Results/Canonical Transformation Dynamical Gutzwiller");
    string loadfunc = "LoadSubResult";
    //        path resdir("/Users/Abuenameh/Documents/Simulation Results/Dynamical Gutzwiller Hartmann Comparison");
#endif
    if (!exists(resdir)) {
        cerr << "Results directory " << resdir << " does not exist!" << endl;
        exit(1);
    }
    ostringstream oss;
    oss << "res." << resi << ".txt";
    path resfile = resdir / oss.str();
    while (exists(resfile)) {
        resi++;
        oss.str("");
        oss << "res." << resi << ".txt";
        resfile = resdir / oss.str();
    }
    args.resi = resi;
    
    cout << "Res: " << resi << endl;

    filesystem::ofstream os(resfile);
    printMath(os, "seed", resi, seed);
    printMath(os, "Delta", resi, D);
    printMath(os, "Wres", resi, Wi);
    printMath(os, "mu0res", resi, mu);

    printMath(os, "tauires", resi, taui);
    printMath(os, "taufres", resi, tauf);
    printMath(os, "ntausres", resi, ntaus);

    printMath(os, "Wires", resi, Wi);
    printMath(os, "Wfres", resi, Wf);

    double m = ntaus / pow(numthreads, 0.5);
    vector<double> tauis, taufs;
    vector<double> ntauss;
    vector<int> tauiis;
    vector<int> taufis;
    tauiis.push_back(0);
    tauis.push_back(taui);
    for (int i = 1; i < numthreads; i++) {
        int j = (int)floor(0.5 + m * pow(i, 0.5));
        if (j >= ntaus)
            break;
        if (j == tauiis.back())
            continue;
        ntauss.push_back(j - tauiis.back());
        taufis.push_back(j-1);
        tauiis.push_back(j);
        tauis.push_back(taui + (tauf - taui)*j/(ntaus-1));
        taufs.push_back(taui + (tauf - taui)*(j-1)/(ntaus-1));
    }
    ntauss.push_back(ntaus - tauiis.back());
    taufis.push_back(ntaus - 1);
    taufs.push_back(tauf);
    
//    cout << tauiis << endl;
//    cout << taufis << endl;
//    cout << ntauss << endl;
//    cout << tauis << endl;
//    cout << taufs << endl;
    
    int len = ntauss.size();
    
    thread_group threads;
    for (int i = 0; i < len; i++) {
        threads.create_thread(bind(&threadfunc, i, tauis[i], taufs[i], ntauss[i], args, resdir));
    os << loadfunc << "[" << resi << "," << i << "];" << endl;
    }
    threads.join_all();

//    os << qwe << "[" << resi << "]=Join@@Table[" << qwe << "[" << resi << ",i],{i,1," << num << "}];";
    combine(os, "taures", resi, len);
    combine(os, "Eires", resi, len);
    combine(os, "Efres", resi, len);
    combine(os, "Qres", resi, len);
    combine(os, "pres", resi, len);
    combine(os, "U0res", resi, len);
    combine(os, "J0res", resi, len);
    combine(os, "b0res", resi, len);
    combine(os, "bfres", resi, len);
    combine(os, "f0res", resi, len);
    combine(os, "ffres", resi, len);
    combine(os, "runtime", resi, len);
//    combine(os, "", resi, len);
    
    os << "subtotalruntime[" << resi << "]=Table[subtotalruntime[" << resi << ",i],{i,0," << len-1 << "}];" << endl;

    ptime end = microsec_clock::local_time();
    time_period period(begin, end);
    cout << endl << period.length() << endl << endl;

    os << "totalruntime[" << resi << "]=\"" << period.length() << "\";" << endl;

    return 0;
}

