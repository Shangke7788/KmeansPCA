#include "tfidf.h"
#include "trie.h"

#include "constant.h"
#include "files.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <map>
#include <math.h>

using namespace std;

int kmp(char* s, char* cmp) {
	static char* last = NULL;
	static int next[1 << 10];
	if (last != cmp) {
		int i = 0, j = -1;
		next[0] = -1;
		while (cmp[i]) {
			if (j == -1 || cmp[i] == cmp[j]) {
				++i, ++j;
				next[i] = j;
			} else {
				j = next[j];
			}
		}
		last = cmp;
	}
	int i = 0, j = 0, res = 0;
	while (s[i]) {
		if (j == -1 || s[i] == cmp[j] || s[i] - 'A' + 'a' == cmp[j]) {
			++i, ++j;
			if (!cmp[j]) {
				res++;
				j = next[j];
			}
		} else {
			j = next[j];
		}
	}
	return res;
}

int word_num(char* s, char* cmp) {
	char *ss = s;
	while (*s) {
		if (*s >= 'A' && *s <= 'Z') {
			*s = *s - 'A' + 'a';
		}
		s++;
	}
	s = ss;
	int res = 0;
	while (*s) {
		if (*s >= 'a' && *s <= 'z') {
			ss = cmp;
			while (*s == *ss) {
				ss++;
				s++;
			}
			if (*s >= 'a' && *s <= 'z') {
				while (*s >= 'a' && *s <= 'z') {
					s++;
				}
				s++;
			} else if (*ss == 0) {
				res++;
				s++;
			} else {
				s++;
			}
		} else {
			s++;
		}
	}
	return res;
}

int document_word_num(const char* filename, char* word) {
	static char s[1 << 10];
	FILE * f;
	f = fopen(filename, "r");
	int num = 0;
	while (fscanf(f, "%s", s) != EOF) {
		//num += word_num(s, word);
		num += kmp(s, word);
	}
	fclose(f);
	return num;
}

void insert_words(Trie* root, const char* filename) {
	static char s[1 << 10];
	FILE * f;
	f = fopen(filename, "r");
	while (fscanf(f, "%s", s) != EOF) {
		char* ss = s;
		while (*ss) {
			if ((*ss >= 'a' && *ss <= 'z') || (*ss >= 'A' && *ss <= 'Z')) {
				ss = root->insert(ss);
			} else {
				ss++;
			}
		}
	}
	fclose(f);
}

Trie* get_words(const char* filename) {
	static char s[1 << 10];
	Trie* r = new Trie();
	FILE *f;
	f = fopen(filename, "r");
	while (fscanf(f, "%s", s) != EOF) {
		char *ss = s;
		while (*ss) {
			if ((*ss >= 'a' && *ss <= 'z') || (*ss >= 'A' && *ss <= 'Z')) {
				if (r->query(ss)) {
					while ((*ss >= 'a' && *ss <= 'z') || (*ss >= 'A' && *ss <= 'Z')) {
						ss++;
					}
				} else {
					ss = r->insert(ss);
				}
			} else {
				ss++;
			}
		}
	}
	fclose(f);
	return r;
}

void get_tfidf_files(const char dic[][50], int num, const char* filename) {
	static char s[1 << 10];
	string f = filename, fout;
	Files fs, tmp;
	fs.clear();
	for (int i = 0; i < num; i++) {
		tmp.getFiles(dic[i]);
		for (Files::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
			fs.push_back(*iter);
		}
	}
	Trie* root = new Trie();

	fprintf(stderr, "Start init... There are 4 steps.\n");
	fprintf(stderr, "step #1: \nRead files starts... This will be take about 1 minute.\n");
	for (int i = 0; i < (int)fs.size(); i++) {
		insert_words(root, fs[i].c_str());
	}
	fprintf(stderr, "Read files ends.\n");
	vector< pair<int, string> > allstring = root->get_all();
	sort(allstring.begin(), allstring.end());
	delete root;

	fout = "tf_" + f + ".txt";
	FILE * tf_all;
	tf_all = fopen(fout.c_str(), "w");
	fprintf(stderr, "step #2: \nWrite the tf.txt file starts...\n");
	for (int i = 0; i < (int)allstring.size(); i++) {
		int y = i;
		strcpy(s, allstring[y].second.c_str());
		fprintf(tf_all, "%-20s %d\n", s, allstring[y].first);
	}
	fclose(tf_all);
	fprintf(stderr, "Write the tf.txt file ends.\n");

	root = new Trie();
	fprintf(stderr, "step #3: \nRead files starts... This will be take about 1 minute.\n");
	for (int i = 0; i < (int)fs.size(); i++) {
		Trie* tmp = get_words(fs[i].c_str());
		*root += *tmp;
		delete tmp;
	}
	fprintf(stderr, "Read files ends.\n");
	allstring = root->get_all();
	sort(allstring.begin(), allstring.end());
	delete root;

	fout = "idf_" + f + ".txt";
	FILE * idf_all;
	idf_all = fopen(fout.c_str(), "w");
	fprintf(stderr, "step #4: \nWrite the idf.txt file starts...\n");
	for (int i = 0; i < (int)allstring.size(); i++) {
		strcpy(s, allstring[i].second.c_str());
		fprintf(idf_all, "%-20s %d\n", s, allstring[i].first);
	}
	fclose(idf_all);
	fprintf(stderr, "Write the idf.txt file ends.\n");

	allstring.clear(), fs.clear();
}

void get_1000_words(const char dic[][50], int snum, const char* filename) {
	static char s[1 << 10];
	string f = filename, fout;
	Files fs, tmp;
	fs.clear();
	for (int i = 0; i < snum; i++) {
		tmp.getFiles(dic[i]);
		for (Files::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
			fs.push_back(*iter);
		}
	}

	FILE *tf_all, *idf_all, *tfidf_all;
	fout = "tf_" + f + ".txt";
	tf_all = fopen(fout.c_str(), "r");
	fout = "idf_" + f + ".txt";
	idf_all = fopen(fout.c_str(), "r");
	if (tf_all == NULL || idf_all == NULL) {
		get_tfidf_files(dic, snum, filename);
		fout = "tf_" + f + ".txt";
		tf_all = fopen(fout.c_str(), "r");
		fout = "idf_" + f + ".txt";
		idf_all = fopen(fout.c_str(), "r");
	}

	vector< pair<double, string> > allstring;
	map<string, int> mp;
	mp.clear(), allstring.clear();

	fprintf(stderr, "Read tf.txt starts...\n");
	int tf_val;
	while (fscanf(tf_all, "%s%d", s, &tf_val) != EOF) {
		string tmp = s;
		mp[tmp] = (int)allstring.size();
		allstring.push_back(make_pair((double)tf_val, tmp));
	}
	fclose(tf_all);
	fprintf(stderr, "Read tf.txt ends.\n");

	fprintf(stderr, "Read idf.txt starts...\n");
	int idf_val;
	while (fscanf(idf_all, "%s%d", s, &idf_val) != EOF) {
		string tmp = s;
		int i = mp[tmp];
		allstring[i].first *= log((double)fs.size() / (double)idf_val);
		mp[tmp] = idf_val;
	}
	fclose(idf_all);
	fprintf(stderr, "Read idf.txt ends.\n");

	fprintf(stderr, "Write tfidf.txt starts...\n");
	sort(allstring.begin(), allstring.end());
	reverse(allstring.begin(), allstring.end());
	fout = "tfidf_" + f + ".txt";
	tfidf_all = fopen(fout.c_str(), "w");
	for (int i = 0; i < (int)allstring.size(); i++) {
		strcpy(s, allstring[i].second.c_str());
		fprintf(tfidf_all, "%-20s %.9lf\n", s, allstring[i].first);
	}
	fclose(tfidf_all);
	fprintf(stderr, "Write tfidf.txt ends...\n");

	fprintf(stderr, "Put 1000 words to Trie and create 1000words.txt...\n");
	fout = "1000words_" + f + ".txt";
	FILE * _1000words;
	_1000words = fopen(fout.c_str(), "w");
	int num = 0;
	for (int i = 0; i < (int)allstring.size(); i++) {
		if (num >= 1000) {
			break;
		}
		strcpy(s, allstring[i].second.c_str());
		if (strlen(s) <= 3) {
			continue;
		}
		if (mp[allstring[i].second] >= 250 * snum || mp[allstring[i].second] < snum * 5) {
			continue;
		}
		fprintf(_1000words, "%-20s %.9lf %d\n", s, log((double)fs.size() / (double)mp[allstring[i].second]), mp[allstring[i].second]);
		num++;
	}
	fclose(_1000words);
	fprintf(stderr, "Now we have 1000 words and 1000words.txt.\n");
}
