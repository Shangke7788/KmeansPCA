#ifndef _SHANGKE_FILES_H_
#define _SHANGKE_FILES_H_

#include <vector>
#include <string>
#include <string.h>
#include <io.h>

using namespace std;

class Files: public vector<string> {
	public:
		/* *
		 * 得到一个目录下的所有文件
		 */
		void getFiles(string path, bool needclear = true) {
			if (needclear) {
				this->clear();
			}
			int hFile = 0;
			struct _finddata_t fileinfo;
			string p;
			if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1) {
				do {
					if (fileinfo.attrib & _A_SUBDIR) {
						if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0) {
							getFiles(p.assign(path).append("\\").append(fileinfo.name), false);
						}
					} else {
						push_back(p.assign(path).append("\\").append(fileinfo.name));
					}
				} while (_findnext(hFile, &fileinfo) == 0);
				_findclose(hFile);
			}
		}
};

const char DICNAME[21][50] = {
	"..\\20_newsgroups",
	"..\\20_newsgroups\\alt.atheism",
	"..\\20_newsgroups\\comp.graphics",
	"..\\20_newsgroups\\comp.os.ms-windows.misc",
	"..\\20_newsgroups\\comp.sys.ibm.pc.hardware",
	"..\\20_newsgroups\\comp.sys.mac.hardware",
	"..\\20_newsgroups\\comp.windows.x",
	"..\\20_newsgroups\\misc.forsale",
	"..\\20_newsgroups\\rec.autos",
	"..\\20_newsgroups\\rec.motorcycles",
	"..\\20_newsgroups\\rec.sport.baseball",
	"..\\20_newsgroups\\rec.sport.hockey",
	"..\\20_newsgroups\\sci.crypt",
	"..\\20_newsgroups\\sci.electronics",
	"..\\20_newsgroups\\sci.med",
	"..\\20_newsgroups\\sci.space",
	"..\\20_newsgroups\\soc.religion.christian",
	"..\\20_newsgroups\\talk.politics.guns",
	"..\\20_newsgroups\\talk.politics.mideast",
	"..\\20_newsgroups\\talk.politics.misc",
	"..\\20_newsgroups\\talk.religion.misc"
};

#endif
