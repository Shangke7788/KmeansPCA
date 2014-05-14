#include "tfidf.h"
#include "trie.h"

#include "constant.h"
#include "files.h"

#include <stdio.h>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <time.h>

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
		if (s[i] == cmp[j]) {
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

int document_word_num(const char* filename, char* word) {
	static char s[1 << 10];
	FILE * f;
	f = fopen(filename, "r");
	int num = 0;
	while (fscanf(f, "%s", s) != EOF) {
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

Trie* get_1000_words() {
	static char s[1 << 10];
	Files fs;
	fs.getFiles(DICNAME[0]);
	Trie* root = new Trie();
	for (int i = 0; i < (int)fs.size(); i++) {
		insert_words(root, fs[i].c_str());
		if ((i + 4) % 100 == 0) {
			printf("%d ok!\n", i);
		}
	}
	vector< pair<int, string> > allstring = root->get_all();
	sort(allstring.begin(), allstring.end());
	delete root;
	root = new Trie();
	int s2 = 110000;
	srand(time(NULL));
	freopen("out.txt", "w", stderr);
	for (int i = 0; i < 1000; i++) {
		int y = rand() % (allstring.size() - s2) + s2;
		strcpy(s, allstring[y].second.c_str());
		if (root->query(s)) {
			i--;
			continue;
		}
		root->insert(s);
		fprintf(stderr, "%d %d %s\n", y, allstring[y].first, s);
	}
	/*
	for (int i = 0; i < 500; i++) {
		strcpy(s, allstring[mid - i].second.c_str());
		root->insert(s);
		fprintf(stderr, "%s %d\n", s, allstring[i].first);
		strcpy(s, allstring[mid + 1 + i].second.c_str());
		root->insert(s);
		fprintf(stderr, "%s %d\n", s, allstring[i].first);
	}
	*/
	return root;
}
