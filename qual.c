#include <stdio.h>
#include <string.h>

#define QUAL_LEN 1000

float phred33ave(char *c) {
	int i, sum = 0, length = strlen(c);
	for (i = 0; i < length; i++) {
		int phred = (int) c[i] - 33;
		sum += phred;
	}
	return sum / length;
}

int main() {

	char qual[QUAL_LEN+1];
	int i, qualsum;
	float avequal;

	printf("seq quality\n");

	scanf("%s", &qual);

	for (i = 0; i < strlen(qual); i++) {
		int phred = (int) qual[i] - 33;
		printf("%c phred score is: %d\n", qual[i], phred);
	}

	float ave = phred33ave(qual);

	printf("average score is %f\n", ave);
	printf("%s\n", qual);

	return 0;
}
