#include<iostream>
#include<cstring>
using namespace std;

std::string exec(const char* cmd) {
	/*
	This function executes an external shell command and retrieves its output in a string.
	Input:
	-cmd : The command to execute.
	Output:
	-result : A string holding the command output.
	*/
	char buffer[128];
	std::string result = "";
#ifdef linux
	FILE* pipe = popen(cmd, "r");
#else
	FILE* pipe = _popen(cmd, "r");
#endif
	if (!pipe) throw std::runtime_error("popen() failed!");
	try {
		while (!feof(pipe)) {
			if (fgets(buffer, 128, pipe) != NULL)
				result += buffer;
		}
	}
	catch (...) {
#ifdef linux
		pclose(pipe);
#else			
		_pclose(pipe);
#endif
		throw;
	}
#ifdef linux
	pclose(pipe);
#else			
	_pclose(pipe);
#endif
	return result;
}

int main(int argc, char* argv[]){
	char cmd[128];
	string n;
#ifdef linux
	strcpy(cmd, "mpirun -n ");
	n=exec("nproc --all");
#else
	strcpy(cmd, "mpiexec -n ");
	n=exec("echo %NUMBER_OF_PROCESSORS%");
#endif
	strcat(cmd, n.c_str());
	for(int i=0; i<strlen(cmd); i++)
		if(cmd[i]=='\r'||cmd[i]=='\n')
			cmd[i]='\0';
	strcat(cmd, " gclust");
	for(int i=1; i<argc; i++){
		strcat(cmd, " ");
		strcat(cmd, argv[i]);
	}
	n=exec(cmd);
	cout<<n.c_str();
	return 0;
}