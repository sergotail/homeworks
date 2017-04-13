# include <stdio.h>
# include <stdlib.h>
# include <inttypes.h>
# include <string.h>

typedef enum { SUCCESS, MEMORY_ERROR, RT_ERROR, EMPTY_ARR, ARG_ERROR } status;

void print_err()
{
	printf("[error]");
}

void free_arr(int64_t ** const arr)
{
	if (*arr) {
		free(*arr);
		*arr = NULL;
	}
	return;
}

size_t min_index(int64_t const * const arr, size_t const arr_size)
{
	size_t res = 0;
	for (size_t i = 1; i < arr_size; ++i)
		if (arr[i] < arr[res])
			res = i;
	return res;
}

void print_arr(int64_t const * const arr, size_t const arr_size)
{
	for (size_t i = 0; i < arr_size; ++i)
		printf("%"PRId64 " ", arr[i]);
	printf("\n");
	return;
}

status N_max(int64_t ** const res, int64_t const * const arr, size_t const arr_size, size_t N)
{
	*res = (int64_t *)malloc((size_t)N * sizeof(int64_t));
	if (!*res)
		return MEMORY_ERROR;
	memcpy(*res, arr, N * sizeof(int64_t));
	if (!*res)
		return MEMORY_ERROR;

	for (size_t i = N; i < arr_size; ++i)
	{
		size_t res_min_index = min_index(*res, N);
		if (arr[i] > (*res)[res_min_index])
			(*res)[res_min_index] = arr[i];
	}

	// sort for decrease order
	for (size_t i = 0; i < N; ++i)
		for (size_t j = N - 1; j > i; j--) {
			if ((*res)[j - 1] < (*res)[j]) {
				int64_t tmp = (*res)[j - 1];
				(*res)[j - 1] = (*res)[j];
				(*res)[j] = tmp;
			}
		}

	return SUCCESS;
}

status read_params(int64_t * const N, int64_t ** const arr, int64_t * const M) {
	if (scanf("%"SCNd64, N) != 1 || *N < 0)
		return ARG_ERROR;
	if (*N)
	{
		if (!((*arr) = (int64_t *)malloc((size_t)(*N) * sizeof(int64_t))))
			return MEMORY_ERROR;

		size_t i = 0;
		while (i < (*N) && scanf("%"SCNd64, &((*arr)[i++])) == 1);
		if (i != *N) {
			free_arr(arr);
			return ARG_ERROR;
		}
	}
	if (scanf("%"SCNd64, M) != 1 || (*N && *M > *N) || (!*N && *M) || (*M < 0)) {
		free_arr(arr);
		return ARG_ERROR;
	}
	if (!*M) {
		if (*N) {
			free_arr(arr);
		}
		return EMPTY_ARR;
	}
	return SUCCESS;
}

int main(int argc, char*argv[])
{
	int64_t N = 0;
	int64_t M = 0;
	int64_t * arr = NULL;
	status st = read_params(&N, &arr, &M);
	if (st != SUCCESS) {
		if (st == EMPTY_ARR)
			return SUCCESS;
		else
			print_err();
			return SUCCESS;
	}
	int64_t * res = NULL;
	if ((st = N_max(&res, arr, (size_t)N, (size_t)M)) != SUCCESS) {
		free_arr(&arr);
		print_err();
		return SUCCESS;
	}
	print_arr(res, (size_t)M);
	free_arr(&res);
	free_arr(&arr);
	return 0;
}
