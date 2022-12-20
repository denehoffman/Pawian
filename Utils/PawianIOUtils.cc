#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "Utils/PawianIOUtils.hh"
#include "ErrLogger/ErrLogger.hh"

namespace PawianIOUtils{
  std::string getFileName(std::string prefix, std::string path){
    std::string  currentCompletePath=prefix+path;
    //check if file exists
    if(!boost::filesystem::exists(currentCompletePath)){
      WarningMsg << "file " << currentCompletePath << " does not exist. Try now to search in absolute path" << endmsg;
      currentCompletePath=path;
      if(!boost::filesystem::exists(currentCompletePath)){
	Alert << "file: " << currentCompletePath << " does not exist!!!!" << endmsg;
	exit(1);
      }
    }
    return currentCompletePath;
  }
}
