// Bioinfo_p1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <new>
#include <stdexcept>
#include <string>
#include <iostream>
#include <queue>
#include <fstream>
#include <sstream>

using namespace std;

//Set up all variables and pointers.
ifstream s1;
ifstream s2;
string ns1;
string ns2;
string seq1;
string seq2;
string labels;
string align1;
string align2;
int** Sub_scores;
double** prob_mat;
int** sub_mat;
int lf;
int sgf;
int Ai;
int Aj;

//Amino acid frequencies.
double AA_freq[] = { 0.08768333990119, 0.0405129182960021, 0.0408784631518651,
		0.0477160345974603, 0.0324709539656211, 0.0378461268859548, 0.0504933695605074,
		0.0898249006830963, 0.0328588505954496, 0.0357514442352249, 0.0852464099207531,
		0.0791031344407513, 0.0148824394639692, 0.0410010190895639, 0.0515802694709073,
		0.0697549720598532, 0.0583275704247605, 0.00931264523877659, 0.0317154088087089,
		0.0630397292098708 };

//Read the amino acid sequences from the input rfiles.
void read_seq_files(string file1, string file2) {
	string bufr;

	s1.open(file1);
	if (s1.fail()) {
		cout << "Failed to open file " << file1 << "\n";
		exit(1);
	}

	s2.open(file2);
	if (s2.fail()) {
		cout << "Failed to open file " << file2 << "\n";
		exit(1);
	}
	//Ignore header
	getline(s1, bufr);
	while (getline(s1, bufr)) {
		seq1.append(bufr);
	}

	// Ignore header
	getline(s2, bufr);
	while (getline(s2, bufr)) {
		seq2.append(bufr);
	}
	//ensure longest sequence is stored as seq_1
	if (seq1.length() > seq2.length()) {
		string temp = seq1;
		seq1 = seq2;
		seq2 = temp;
	}

	//Create a 2D array for substitution scores.
	Sub_scores = new int*[ seq1.length() + 1];
	for (int i = 0; i < (seq1.length() + 1); i++) {
		Sub_scores[i] = new int[seq2.length() + 1];
	}
}

//Find the substitution score.
int find_s_val(int i, int j) {
	int si1;
	int si2;
	for (int h = 0; h < 20; h++) {
		if (seq1[i] == labels[h]) {
			si1 = h;
			break;
		}
	}
	for (int k = 0; k < 20; k++) {
		if (seq2[j] == labels[k]) {
			si2 = k;
			break;
		}
	}
	return (sub_mat[si1][si2]);
}

void pam_read(string if3, int divergence) {
	double** temp_mat = new double*[20];
	for (int l = 0; l < 20; l++) {
		temp_mat[l] = new double[20];
	}
	double** pow_mat = new double* [20];
	for (int l = 0; l < 20; l++) {
		pow_mat[l] = new double[20];
	}
	prob_mat = new double* [20];
	for (int l = 0; l < 20; l++) {
		prob_mat[l] = new double[20];
	}
	string line_buf;
	ifstream pf;
	pf.open(if3);
	if (pf.fail()) {
		cout << "Failed to open file " << if3 << "\n";
		exit(1);
	}
	//ignore header
	getline(pf, line_buf);
	getline(pf, line_buf);
	//tokenize line_buf for labels
	stringstream tmp_stream(line_buf);
	string tmp;
	for (int lab_ind = 0; lab_ind < 20; lab_ind++) {
		getline(tmp_stream, tmp, ',');
		labels.push_back(tmp[0]);
	}
	//int column_total [20];
	//get first line of probabilities and tokenize each value from there
	for (int i = 0; i < 20; i++) {
		getline(pf, line_buf);
		//column_total[i] = 0;
		string t;
		stringstream t_stream(line_buf);
		for (int j = 0; j < 20; j++) {
			getline(t_stream, t, ',');
			temp_mat[i][j] = ((double) atof(t.c_str()) / (double) 10000);
			prob_mat[i][j] = temp_mat[i][j];
		}
	}


	for (int d = 2; d <= divergence; d++) {
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++) {
				double sum = 0;
				for (int z = 0; z < 20; z++) {
					sum = sum + (prob_mat[i][z] * temp_mat[z][j]);
				}
				pow_mat[i][j] = sum;
			}
		}
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++) {
				prob_mat[i][j] = pow_mat[i][j];
			}
		}
	}
	cout << "Frequencies:\n";
	for (int i = 0; i < 20; i++) {
		cout << labels[i] << " ";
	}
	cout << "\n";
	for (int i = 0; i < 20; i++) {
		cout << labels[i] << ": ";
		for (int j = 0; j < 20; j++) {
			cout << prob_mat[i][j] << " ";
		}
		cout << "\n";
	}
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			sub_mat[i][j] = sub_mat[j][i] = (int)(0.5 + (10 * (log(prob_mat[j][i] / AA_freq[j])) / log(10) + 10 * (log(prob_mat[i][j] / AA_freq[i])) / log(10))/2);
		}
	}
	cout << "substitution matrix:\n";
	for (int i = 0; i < 20; i++) {
		cout << labels[i] << " ";
	}
	cout << "\n";
	for (int i = 0; i < 20; i++) {
		cout << labels[i] << ": ";
		for (int j = 0; j < 20; j++) {
			cout << sub_mat[i][j] << " ";
		}
		cout << "\n";
	}
}

void sub_read(string if3) {
	string line_buf;
	ifstream pf;
	pf.open(if3);
	if (pf.fail()) {
		cout << "Failed to open file " << if3 << "\n";
		exit(1);
	}
	//ignore header
	getline(pf, line_buf);
	getline(pf, line_buf);
	//tokenize line_buf for labels
	stringstream tmp_stream(line_buf);
	string tmp;
	for (int lab_ind = 0; lab_ind < 20; lab_ind++) {
		getline(tmp_stream, tmp, ',');
		labels.push_back(tmp[0]);
	}
	for (int i = 0; i < 20; i++) {
		getline(pf, line_buf);
		string t;
		stringstream t_stream(line_buf);
		for (int j = 0; j < 20; j++) {
			getline(t_stream, t, ',');
			sub_mat[i][j] = atoi(t.c_str());
		}
	}

}

void Global(int g) {
	for (int i = 0; i < (seq1.length() + 1); i++) {
		for (int j = 0; j < (seq2.length() + 1); j++) {
			Sub_scores[i][j] = 0;
		}
	}
	for (int i = 0; i < (seq1.length() + 1); i++) {
		Sub_scores[i][0] = i * g;
	}
	for (int j = 0; j < (seq2.length() + 1); j++) {
		Sub_scores[0][j] = j * g;
	}
	for (int i = 1; i < (seq1.length() + 1); i++) {
		for (int j = 1; j < (seq2.length() + 1); j++) {
			int s = find_s_val((i-1), (j-1));
			Sub_scores[i][j] = max({ Sub_scores[i - 1][j] + g,
				Sub_scores[i][j - 1] + g, Sub_scores[i - 1][j - 1] + s });
		}
	}
}

void Local(int g) {
	for (int i = 0; i < (seq1.length() + 1); i++) {
		for (int j = 0; j < (seq2.length() + 1); j++) {
			Sub_scores[i][j] = 0;
		}
	}
	for (int i = 0; i < (seq1.length() + 1); i++) {
		Sub_scores[i][0] = 0;
	}
	for (int j = 0; j < (seq2.length() + 1); j++) {
		Sub_scores[0][j] = 0;
	}
	for (int i = 1; i < (seq1.length() + 1); i++) {
		for (int j = 1; j < (seq2.length() + 1); j++) {
			int s = find_s_val((i-1), (j - 1));
			Sub_scores[i][j] = max({ Sub_scores[i - 1][j] + g,
				Sub_scores[i][j - 1] + g, Sub_scores[i - 1][j - 1] + s, 0});
		}
	}
}

void Semi_Global(int g) {
	for (int i = 0; i < (seq1.length() + 1); i++) {
		for (int j = 0; j < (seq2.length() + 1); j++) {
			Sub_scores[i][j] = 0;
		}
	}

	for (int i = 1; i < (seq1.length() + 1); i++) {
		for (int j = 1; j < (seq2.length() + 1); j++) {
			int s = find_s_val((i - 1), (j - 1));
			Sub_scores[i][j] = max({ Sub_scores[i - 1][j] + g,
				Sub_scores[i][j - 1] + g, Sub_scores[i - 1][j - 1] + s });
		}
	}
}

void Align(int i, int j, int g) {
	string st;
	if ((i == 0) && (j == 0)) {
		align1.insert(0, "_");
		align2.insert(0, "_");
		return;
	}
	else {
		int v = Sub_scores[i - 1][j] + g;
		int h = Sub_scores[i][j - 1] + g;
		if (Sub_scores[i][j] == v) {
			if ((Sub_scores[i][j] != 0) && (lf == 1)) {
				align2.insert(0, "_");
				align1.insert(0, 1, seq1[i-1]);
				Align(i - 1, j, g);
			}
			else if (lf == 0) {
				align2.insert(0, "_");
				align1.insert(0, 1, seq1[i-1]);
				Align(i - 1, j, g);
			}
			else {
				align2.insert(0, "_");
				align1.insert(0, 1, seq1[i-1]);
			}
			
		}
		else if(Sub_scores[i][j] == h){
			if ((Sub_scores[i][j] != 0) && (lf == 1)) {
				char t = seq2[j-1];
				align1.insert(0, "_");
				align2.insert(0, 1, t);
				Align(i, j - 1, g);
			}
			else if (lf == 0) {
				char t = seq2[j-1];
				align1.insert(0, "_");
				align2.insert(0, 1, t);
				Align(i, j - 1, g);
			}
			else {
				char t = seq2[j-1];
				align1.insert(0, "_");
				align2.insert(0, 1, t);
			}
			
		}
		else {
			if ((Sub_scores[i][j] != 0) && (lf == 1)) {
				align1.insert(0, 1, seq1[i-1]);
				align2.insert(0, 1, seq2[j-1]);
				Align(i-1, j - 1, g);
			}
			else if (lf == 0) {
				align1.insert(0, 1, seq1[i-1]);
				align2.insert(0, 1, seq2[j-1]);
				Align(i-1, j - 1, g);
			}
			else {
				align1.insert(0, 1, seq1[i-1]);
				align2.insert(0, 1, seq2[j-1]);
			}
			
		}
	}
}
//find where to start alignment
void find_align_in() {
	//search for max in last column and row for semi-global
	if (sgf == 1) {
		int max_c = INT_MIN;
		int max_r = INT_MIN;
		for (int i = 0; i < (seq1.length() + 1); i++) {
			if (max_c < Sub_scores[i][seq2.length()]) {
				max_c = Sub_scores[i][seq2.length()];
				Ai = i;
			}
		}
		for (int j = 0; j < (seq2.length() + 1); j++) {
			if (max_r < Sub_scores[seq1.length()][j]) {
				max_r = Sub_scores[seq1.length()][j];
				Aj = j;
			}
		}
		if (max_c > max_r) {
			Aj = seq2.length();
		}
		else {
			Ai = seq1.length();
		}
	}
	//find max of all values for local
	else if (lf == 1) {
		int max_val = INT_MIN;
		for (int i = 0; i < (seq1.length() + 1); i++) {
			for (int j = 0; j < (seq2.length() + 1); j++) {
				if (Sub_scores[i][j] > max_val) {
					max_val = Sub_scores[i][j];
					Ai = i;
					Aj = j;
				}
			}
		}
	}
	else {
		Ai = seq1.length();
		Aj = seq2.length();
	}
}

string translate(string codon) {
	string amino;
	if (codon == "GCU" || "GCC" || "GCA" || "GCG" || "GCT") {
		amino = "A";
	}
	else if (codon == "CGU" || "CGC" || "CGA" || "CGG" || "CGT" || "AGA" || "AGG") {
		amino = "R";
	}
	else if (codon == "AAU" || "AAC" || "AAT") {
		amino = "N";
	}
	else if (codon == "GAU" || "GAC" || "GAT") {
		amino = "D";
	}
	else if (codon == "UGU" || "UGC" || "TGT" || "TGC") {
		amino = "C";
	}
	else if (codon == "CAA" || "CAG") {
		amino = "Q";
	}
	else if (codon == "GAA" || "GAG") {
		amino = "E";
	}
	else if (codon == "GGU" || "GGC" || "GGA" || "GGG" || "GGT") {
		amino = "G";
	}
	else if (codon == "CAU" || "CAC" || "CAT") {
		amino = "H";
	}
	else if (codon == "AUU" || "AUC" || "AUA" || "ATT" || "ATC" || "ATA") {
		amino = "I";
	}
	else if (codon == "UUA" || "UUG" || "TTA" || "TTG" || "CUU" || "CUA" || "CUC" || "CUG" || "CTT" || "CTA" || "CTC" || "CTG") {
		amino = "L";
	}
	else if (codon == "AAA" || "AAG") {
		amino = "K";
	}
	else if (codon == "AUG" || "ATG") {
		amino = "M";
	}
	else if (codon == "UUU" || "UUC" || "TTT" || "TTC") {
		amino = "F";
	}
	else if (codon == "CCU" || "CCC" || "CCA" || "CCG" || "CCT") {
		amino = "P";
	}
	else if (codon == "UCU" || "UCC" || "UCA" || "UCG" || "TCT" || "TCC" || "TCA" || "TCG" || "AGU" || "AGC" || "AGT") {
		amino = "S";
	}
	else if (codon == "ACU" || "ACC" || "ACA" || "ACG" || "ACT") {
		amino = "T";
	}
	else if (codon == "UGG" || "TGG") {
		amino = "W";
	}
	else if (codon == "UAU" || "UAC" || "TAT" || "TAC") {
		amino = "Y";
	}
	else if (codon == "GUU" || "GUC" || "GUA" || "GUG" || "GTT" || "GTC" || "GTA" || "GTG") {
		amino = "V";
	}
	else {
		amino = "X";
	}
	return amino;
}

void read_n(string if1, string if2) {
	string buf;
	s1.open(if1);
	if (s1.fail()) {
		cout << "Failed to open file " << if1 << "\n";
		exit(1);
	}

	s2.open(if2);
	if (s2.fail()) {
		cout << "Failed to open file " << if2 << "\n";
		exit(1);
	}

	getline(s1, buf);
	while (getline(s1, buf)) {
		ns1.append(buf);
	}

	getline(s2, buf);
	while (getline(s2, buf)) {
		ns2.append(buf);
	}

	//ensure longest sequence is row length
	if (ns1.length() > ns2.length()) {
		string temp = ns1;
		ns1 = ns2;
		ns2 = temp;
	}

	for (int i = 0; i < (ns1.length() - 2); i = i + 3) {
		string codon;
		string amino;
		codon.push_back(ns1[i]);
		codon.push_back(ns1[i] + 1);
		codon.push_back(ns1[i] + 2);
		amino = translate(codon);
		seq1.append(amino);
	}
	for (int i = 0; i < (ns2.length() - 2); i = i + 3) {
		string codon;
		string amino;
		codon.push_back(ns2[i]);
		codon.push_back(ns2[i] + 1);
		codon.push_back(ns2[i] + 2);
		amino = translate(codon);
		seq2.append(amino);
	}

	cout << "Translations:\n";
	cout << ns1 << " = " << seq1 << "\n";
	cout << ns2 << " = " << seq2 << "\n";
}

int main()
{
	lf = 0;
	sgf = 0;
	string if1;
	string if2;
	string if3;
	int divergence;
	int g;
	int type;
	char pam_ans;
	char n_ans;
	sub_mat = new int* [20];
	for (int l = 0; l < 20; l++) {
		sub_mat[l] = new int[20];
	}
	
	cout << "Enter the name of the file containing the first sequence.\n";
	cin >> if1;
	cout << "Enter the name of the file containing the sequence to be compared\n";
	cin >> if2;
	cout << "Is this a DNA/RNA sequence? (y or n)\n";
	cin >> n_ans;
	if (n_ans == 'y') {
		read_n(if1, if2);
	}
	else {
		read_seq_files(if1, if2);
	}
	cout << "Enter the name of the file containing a probability or substitution matrix\n";
	cin >> if3;
	cout << "is this a PAM probability matrix? (y or n)\n";
	cin >> pam_ans;
	if (pam_ans == 'y') {
		cout << "Enter the degree of evolutionary divergence\n";
		cin >> divergence;
		pam_read(if3, divergence);
	}
	else {
		sub_read(if3);
	}
	cout << "Enter the number corresponding to the alignment type\n1: global\n2: semi global\n3: local\n";
	cin >> type;
	cout << "Enter the gap score\n";
	cin >> g;
	switch (type) {
	case 1:
		Global(g);
		break;
	case 2:
		sgf = 1;
		Semi_Global(g);
		break;
	case 3:
		lf = 1;
		Local(g);
		break;
	}

	cout << "Substitution scores:\n";
	for (int i = 0; i < (seq2.length() + 1); i++) {
		if (i == 0) {
			cout << "   _: ";
		}
		else {
			cout << seq2[i - 1] << " ";
		}
	}
	cout << "\n";
	for (int i = 0; i < (seq1.length() + 1); i++) {
		if (i == 0) {
			cout << "_: ";
		}
		else {
			cout << seq1[i - 1] << ": ";
		}
		
		for (int j = 0; j < (seq2.length() + 1); j++) {
			cout << Sub_scores[i][j] << " ";
		}
		cout << "\n";
	}
	find_align_in();
	Align(Ai, Aj, g);
	cout << "score: " << Sub_scores[Ai][Aj] << "\n";
	cout << "Alignment:\n";
	cout << align1 << "\n" << align2 << "\n";
}

