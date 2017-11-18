//
//  create-folders.h
//  PBRTV3-EGSR2017
//
//  Created by Gurpreet Singh Bagga on 26/02/17.
//
//

#ifndef PBRTV3_EGSR2017_create_folders_h
#define PBRTV3_EGSR2017_create_folders_h

#include <unistd.h> //to use getcwd(cwd, sizeof(cwd));
#include <sys/stat.h> // mkdir
#include <sstream>


void create_folders(std::string homedir, std::string &data, std::string &images, std::string &graphs){
    
    std::time_t rawtime;
    std::tm* timeinfo;
    char buffer [80];
    
    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    
    std::strftime(buffer,80,"%Y-%m-%d",timeinfo);
    std::puts(buffer);
    
    //#####################Setting up folders to store data#########################
    
    std::stringstream ss, ssdate, ssdir, ssdata, ssimg, ssgraph;
    //ssdate.imbue(std::locale(ssdate.getloc(), facet));
    ssdate << buffer;
    
    //const char* homeDir = getenv ("Home");
    //char final [256];
    //sprintf (final, "%s/Desktop/%s",homeDir, game_name);
    
    ssdir << homedir << ssdate.str();
    if(!mkdir(ssdir.str().c_str(),0775))
        std::cerr <<"Directory with date!"<< std::endl;
    
    ssdata << ssdir.str() << "/datafiles/";
    if(!mkdir(ssdata.str().c_str(),0775))
        std::cerr <<"Directory with datafiles!"<< std::endl;
    
    ssimg << ssdir.str() << "/images/";
    if(!mkdir(ssimg.str().c_str(),0775))
        std::cerr <<"Directory with images!"<< std::endl;
    
    ssgraph << ssdir.str() << "/graphs/";
    if(!mkdir(ssgraph.str().c_str(),0775))
        std::cerr <<"Directory with graphs!"<< std::endl;
    
    data = ssdata.str();
    images = ssimg.str();
    graphs = ssgraph.str();
    
}

#endif
