#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cstring>
#include<sstream>
#include<iomanip>
#include<cmath>
using namespace std;

char progPath[255]=""; //Path of the program directory
string mDist;
string *nc; //Array that will hold the input genomes' names
double **MatSimil, **MD; //MatSimil will hold the similarity matrix and MD may hold a distance matrix input as a txt file
int nbGenomes; //The number of genomes in the input fasta file initialized in similarity() function

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

double ** matriceDistances(string **dicoMuscle, int n = nbGenomes) {
	/*
	This function is a subfunction of "similarity".
	Input:
	-dicoMuscle : A dictionary of sequences in which the keys are the names of the sequences.
	-nc : The names of the sequences.
	Output:
	-Matrix : A matrix of distance.
	*/
	double **matDist;
	matDist = new double*[n];
	for (int i = 0; i < n; i++)
		matDist[i] = new double[n];
	
	for (int i = 0; i<n; i++)
		for (int j = i; j<n; j++){
			char cmd[5000] = "";
			strcpy(cmd,progPath);
			strcat(cmd,"/gdist -mdist ");
			strcat(cmd,mDist.c_str());
			strcat(cmd," -ali1 ");
			strcat(cmd,dicoMuscle[i][1].c_str());
			strcat(cmd," -ali2 ");
			strcat(cmd,dicoMuscle[j][1].c_str());
			matDist[i][j] = stod(exec(cmd));
		}

	for (int i = 0; i<n; i++) 
		for (int j = 0; j<i; j++)
			matDist[i][j] = matDist[j][i];

	return matDist;
}

void listGenomes(string fListe){
	ifstream inFile(fListe.c_str());
	if(!inFile)
	{
		cout<<"Couldn't open input fasta file"<<endl;
		exit(1);
	}
	string line;
	int i=0;
	while(getline(inFile, line)){
		if(line=="")
			continue;
        if(line.at(0)=='>'||line.at(0)==';'){
			line.erase(0,1);
			#ifdef linux
				if(line[strlen(line.c_str())-1]=='\r')
					line[strlen(line.c_str())-1]='\0';
			#endif
			nc[i]=line.c_str();
			i++;
		}
	}
}

void similarity(string fListe) {
	/*
	This function returns the similarity matrix.
	Input:
	-fListe : Name of the fasta file holding a list of sequences with their names.
	Output:
	-MatSimil : The similarity matrix (global variable).
	-nc : Sequential number with the names of the sequences (global variable).
	*/

	// I) First step : Apply

	char cmd[150] = "";
	strcpy(cmd,progPath);
	strcat(cmd,"/muscle -quiet -in ");
	strcat(cmd, fListe.c_str());
	string liste = exec(cmd);
	// Save the result
	/*std::ofstream out("aligned.txt");
	out << liste;
	out.close();*/
	// Read the alignmed result from a saved file
	/*std::ifstream ifs("aligned.txt");
	std::string liste( (std::istreambuf_iterator<char>(ifs) ),
                       (std::istreambuf_iterator<char>()    ) );*/
	int index;
	string **dicoMuscle,tempGen;
	double **MatDistance, max = 0.0;
	nbGenomes=0;
	for (int i = 0; i<liste.size(); i++)
		if ((liste[i] == '>' && liste[i+1] != '>') || (liste[i] == ';' && liste[i+1] != ';')) //Taking into consideration the error of duplication by checking liste[i+1]
			nbGenomes++;
	dicoMuscle = new string*[nbGenomes];
	nc = new string[nbGenomes];
	listGenomes(fListe);
	MatSimil = new double*[nbGenomes];
	for (int i = 0; i<nbGenomes; i++) {
		dicoMuscle[i] = new string[2];
		MatSimil[i] = new double[nbGenomes];
	}
	for (int i = 0; i<liste.size();) {
		if (liste[i] == '>'||liste[i]==';') {
			tempGen="";
			index=0;
			i++;
			do {
				tempGen += liste[i];
				i++;
			} while (liste[i] != '\n'&&liste[i]!='\r');
			while(tempGen.compare(nc[index])!=0){//strcmp(nc[index].c_str(),tempGen.c_str())!=0){
				index++;
			}
			dicoMuscle[index][0]=tempGen;
		}
		else {
			do {
				if (liste[i] != '\n') {
					dicoMuscle[index][1] += liste[i];
				}
				i++;
			} while (liste[i] != '>'&&liste[i] != ';'&&i<liste.size());
		}
	}
	//MatDistance=MD; //Set MatDisctance to the distance matrix MD imported from the external text file
	MatDistance = matriceDistances(dicoMuscle); //Calculates the distance matrix
	
	for (int i = 0; i<nbGenomes; i++)
		for (int j = 0; j<nbGenomes; j++)
			if (MatDistance[i][j]>max)
				max = MatDistance[i][j];
	for (int i = 0; i<nbGenomes; i++)
		for (int j = 0; j<nbGenomes; j++)
			MatDistance[i][j] /= max;
	for (int i = 0; i<nbGenomes; i++)
		for (int j = 0; j<nbGenomes; j++)
			MatSimil[i][j] = 1 - MatDistance[i][j];
	/*	// Display the matrix
	for(int i=0;i<100;i++){
	for(int j=0;j<100;j++){
	std::cout<<setw(16);
	std::cout << std::fixed;
	std::cout << std::setprecision(12);
	std::cout <<MatSimil[i][j];
	}
	std::cout<<endl;
	}*/

	//Free memory from dynamically allocated matrixes
	for (int i = 0; i<nbGenomes; i++) {
		if (MatDistance != MD)
			delete[] MatDistance[i];
		delete[] dicoMuscle[i];
	}
	if (MatDistance != MD)
		delete[] MatDistance;
	delete[] dicoMuscle;
}

string getStringFromDouble(double num){
	ostringstream streamObj;
	streamObj << std::fixed;
	streamObj << std::setprecision(18);
	streamObj << num;
	string strObj = streamObj.str();
	return strObj;
}

int main(int argc, char* argv[]) {
	string fListe,fGroupes;
	string groupes;

	//Displaying usage message in case called without arguments
	if(argc==1)
		cout<<"You launched the program without arguments!\n\n"
		<<"By default this program reads the genomes from the file named 'sequences.fasta' and outputs the clustering in the file 'Clustering.txt' in the current path.\n"
		<<"The default distance matrix used is EDNAFULL.\n\n"
		<<"You can customize this configuration by setting the arguments -in for selecting the input fasta file, -out for selecting the output txt file, and -mdist for selecting a different matrix.\n"
		<<"The currently configured matrices are: EDNAFULL, BLOSUM62, and PAM250."<<endl<<endl;

	//Setting the path of our executable
#ifdef linux
	realpath(argv[0], progPath);
#else			
	_fullpath(progPath, argv[0], sizeof(progPath));
#endif
	
	progPath[strlen(progPath)-7]='\0';
	
	//Setting the default values for the arguments
	fListe = string(progPath);
	fListe += "/sequences.fasta";
	fGroupes = string(progPath);
	fGroupes += "/Clustering.txt";
	mDist="EDNAFULL";

	//Reading the input arguments
	for(int i=1; i<argc; i+=2)
		if(strcmp(argv[i],"-mdist")==0)
			mDist=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(argv[i],"-in")==0)
			fListe=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(argv[i],"-out")==0)
			fGroupes=argv[i+1];

	//Validate the input arguments
	if(mDist!="EDNAFULL" && mDist!="BLOSUM62" && mDist!="PAM250"){
		cout<<"Invalid distance matrix. Available matrices are: EDNAFULL, BLOSUM62, and PAM250.\n";
		return 1;
	}

	if (fListe.substr(fListe.find_last_of(".") + 1) != "fasta" && fListe.substr(fListe.find_last_of(".") + 1) != "dat"){
			cout << "Invalid input filename! Please input a fasta file containing the genomes using the argument -in.\n";
			return 1;
		}

	if (fGroupes.substr(fGroupes.find_last_of(".") + 1) != "txt" && fGroupes.substr(fGroupes.find_last_of(".") + 1) != "dat"){
		cout << "Invalid output filename! Please name a txt file for the results using the argument -out.\n";
		return 1;
	}
	
	//Check if the input fasta file is accessible
	std::ifstream infile(fListe);
    if(!infile.good()){
		cout<<"The input fasta file is missing or not accessible.\n\n";
		return 1;
	}
	
	//Calculate the similarity matrix
	similarity(fListe);
	string matriceSimilitude="";
	string refsGenomes="";
	
	for(int i=0; i<nbGenomes; i++){
		for(int j=0; j<nbGenomes; j++){
			matriceSimilitude+=getStringFromDouble(MatSimil[i][j]);
			matriceSimilitude+=" ";
		}
		matriceSimilitude+="\r\n";
		
		refsGenomes+=nc[i];
		refsGenomes+=",";
	}
#ifdef linux
	//refsGenomes.erase( std::remove(refsGenomes.begin(), refsGenomes.end(), '\0'), refsGenomes.end() );
	//refsGenomes[3*nbGenomes+refsLen]=']';
#else
	//refsGenomes[strlen(refsGenomes.c_str())-1]=']';

#endif		
		
	std::ofstream out1("matSimil.txt");
	out1 << matriceSimilitude;
	out1.close();
	std::ofstream out2("refs.txt");
	out2 << refsGenomes;
	out2.close();
	
	//Free memory from dynamically allocated matrixes	
	for (int i = 0; i<nbGenomes; i++) {
		delete[] MatSimil[i];
		//delete [] MD[i];
	}
	delete[] MatSimil, delete[] MD, delete[] nc;
	
	char cmd[150]="";
	strcat(cmd,progPath);
#ifdef linux
	strcat(cmd,"/gclust-GMM");
#else
	strcpy(cmd,"python ");
	strcat(cmd,progPath);
	strcat(cmd,"/gclust-GMM.py");
#endif
	exec(cmd);
	//remove("matSimil.txt");
	//remove("refs.txt");
	rename("Clustering.txt",fGroupes.c_str());

	return 0;
}