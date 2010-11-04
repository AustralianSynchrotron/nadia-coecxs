#include <iostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "Config.h"

#define FAILURE 0
#define SUCCESS 1

#define KEY 0
#define VALUE 1

using namespace std;

Config::Config(string file_name){
  //open the file and
  //parse it's contents.
     
    ifstream file(file_name.c_str());
    if(!file.is_open()){
      cout << "Could not open the file " << file_name << endl;
      status = FAILURE;
      return;
    }

    string line;
    
    while(!file.eof()){
      getline(file,line);
      
      //tokenise the line of text

      //first remove the comments
      size_t pos = line.find("#");
      if(pos!=string::npos){
	line = line.substr(0,pos);
      }

      //now we get the key part and value part of the string
      pos = line.find("=");
      if(pos!=string::npos){
      string key_temp = line.substr(0,pos);

      //now get rid of white space      
      istringstream key_stream(key_temp);
      string key;      
      key_stream >> key;

      string value = line.substr(pos+1);

      mapping.insert(pair<string,string>(key,value));

      /**
      istringstream value_stream(value_temp);
      string value;
      value_stream >> value;
      
      cout << "original line after comment removed:"<<line<<endl;
      cout << "key with shite space:"<<key_temp<<"-"<<endl;
      cout << "key without shite space:"<<key<<"-"<<endl;
      cout << "value with shite space:"<<value_temp<<"-"<<endl;
      cout << "value without shite space:"<<value<<"-"<<endl; **/
      }
    }
    
    status=SUCCESS;
   
};
 


string Config::getString(string key){
  MapType::iterator iter = mapping.find(key);
  if (iter == mapping.end() ){
    std::cout << "Key is not in my_map" << endl;
    status=FAILURE;
    return "";
  }
  string value;      
  istringstream value_temp(iter->second);
  value_temp >> value;
  //remove white space
  
  return value;
  
}

double Config::getDouble(string key){
MapType::iterator iter = mapping.find(key);
 if (iter == mapping.end() ){
   std::cout << "Key is not in my_map" << endl;
   status=FAILURE;
   return 0;
 }
 return atof((iter->second).c_str());
 
}

int Config::getInt(string key){
MapType::iterator iter = mapping.find(key);
 if (iter == mapping.end() ){
   std::cout << "Key is not in my_map" << endl;
   status=FAILURE;
   return 0;
 }
 return atoi((iter->second).c_str());
 
}

  
