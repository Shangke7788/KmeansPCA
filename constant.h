#ifndef _SHANGKE_CONSTANT_H_
#define _SHANGKE_CONSTANT_H_

#ifdef __cplusplus
extern "C" {
#endif

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

const char A2[2][50] = {
	"..\\20_newsgroups\\alt.atheism",
	"..\\20_newsgroups\\comp.graphics"
};

const char B2[2][50] = {
	"..\\20_newsgroups\\talk.politics.mideast",
	"..\\20_newsgroups\\talk.politics.misc"
};

const char A5[5][50] = {
	"..\\20_newsgroups\\comp.graphics",
	"..\\20_newsgroups\\rec.motorcycles",
	"..\\20_newsgroups\\rec.sport.baseball",
	"..\\20_newsgroups\\sci.space",
	"..\\20_newsgroups\\talk.politics.mideast"
};

const char B5[5][50] = {
	"..\\20_newsgroups\\comp.graphics",
	"..\\20_newsgroups\\comp.os.ms-windows.misc",
	"..\\20_newsgroups\\rec.autos",
	"..\\20_newsgroups\\sci.electronics",
	"..\\20_newsgroups\\talk.politics.misc"
};

const int BALANCE2[2] = { 100, 100 };
const int BALANCE5[5] = { 100, 100, 100, 100, 100 };
const int UNBALANCE5[5] = { 200, 140, 120, 100, 60 };

#ifdef __cplusplus
}
#endif

#endif
