#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <map>
#include <list>

#define FAILURE 0
#define SUCCESS 1

typedef std::map<std::string,std::string> MapType;

class Config{

 private:
  MapType mapping;
  int status;

 public:

  Config(std::string file_name);
  //~Config();

  std::string getString(std::string key);
  double getDouble(std::string key);
  int getInt(std::string key);
  //bool getBoolean(std::string key);

  std::list<int> * getIntList(std::string key);
  std::list<std::string> * getStringList(std::string key);

  int getStatus(){return status;};


 private:
    
};

#endif
