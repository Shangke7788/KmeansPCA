#ifndef _SHANGKE_TFIDF_H_
#define _SHANGKE_TFIDF_H_

#include "trie.h"

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

int kmp(char* s, char* cmp);

int word_num(char* s, char* cmp);

int document_word_num(const char* filename, char* word);

void insert_words(Trie* root, const char* filename);

Trie* get_words(const char* filename);

void get_tfidf_files(const char dic[][50], int n, const char* filename);

void get_1000_words(const char dic[][50], int n, const char* filename);

#ifdef __cplusplus
}
#endif

#endif
