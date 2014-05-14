#ifndef _SHANGKE_TFIDF_H_
#define _SHANGKE_TFIDF_H_

#include "trie.h"

using namespace std;

int kmp(char* s, char* cmp);

int document_word_num(const char* filename, char* word);

void insert_words(Trie* root, const char* filename);

Trie* get_1000_words();

#endif
