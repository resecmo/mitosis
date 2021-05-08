#include "mt.hpp"
#include "bulk.hpp"

using namespace std;


void mt2bulk();

void bulk2mt();

void step_full(){
	step_bulk();//osetr
	mt2bulk();//?order
	step_mt();//basilius
	bulk2mt();//order?
	if (t % each == 0){
		stamp_bulk(outStreamBulk);
		stamp_mt(outStreamMT);
	}
}

void exchange();

void step_full(){
	step_bulk();
	step_mt();
	exchange();
	if (t % each == 0){
		stamp_bulk(outStreamBulk);
		stamp_mt(outStreamMT);
	}
}


int main(){
	//read args
	
	initialize_mt();//basilius
	initialize_bulk();//osetr

	for (long long t = 0; t < T; ++t){
		step_full();
	}


	return 0;
}
