# include <stdlib.h>
# include <stdio.h>
# include <inttypes.h>
# include <ctype.h>
# include <wctype.h>
# include <stdbool.h>
# include <string.h>

typedef enum status { SUCCESS, LEX_ERROR, SYNT_ERROR, RT_ERROR, MEMORY_ERROR } status;
typedef enum lex_type { EMPTY_LEX, EOF_LEX, OPEN, CLOSE, PLUS, MINUS, MUL, DIV, NUMBER } lex_type;
typedef int32_t bigint_digit;
typedef struct lex
{
	lex_type type;
	char * val;
	size_t size;
} lex;
typedef struct bigint
{
	size_t digits_size;
	bigint_digit * digits;
	bool is_negative;
} bigint;

// declarations for long arithmetic
static int8_t const BIGINT_BASE_DEC = 4;
static int16_t const BIGINT_BASE = 10000;

// functions
// common
void print_err();

// for parser
lex create_lex();
void free_lex(lex * const);
status read_lex(lex * const);
void print_lex(lex const * const);
status calc_expression(lex * const, bigint * const);
status expression(lex * const, bigint * const);
status term(lex * const, bigint * const);
status factor(lex * const, bigint * const);

// for long arithmetic
bigint create_bigint();
void free_bigint(bigint * const);
void print_bigint(bigint const * const);
bool is_bigint_zero(bigint const * const);
bool is_left_bigint_ge_abs(bigint const * const, bigint const * const);
status n_digits_from_string(bigint_digit * const, char const * const, size_t const);
status bigint_from_string(bigint * const, char const * const, size_t const);
status bigint_copy(bigint * const, bigint const * const);
status bigint_abs(bigint * const, bigint const * const);
status bigint_delete_lead_zeros(bigint * const);
status bigint_insert(bigint * const, size_t const, bigint_digit const);
status bigint_push_back(bigint * const, bigint_digit const);
status bigint_push_front(bigint * const, bigint_digit const);
status bigint_nonzero_nonneg_plus(bigint * const, bigint const * const, bigint const * const);
status bigint_nonzero_nonneg_minus(bigint * const, bigint const * const, bigint const * const);
status bigint_nonzero_nonneg_mul(bigint * const, bigint const * const, bigint const * const);
status bigint_nonzero_nonneg_div(bigint * const, bigint const * const, bigint const * const);
status bigint_plus(bigint * const, bigint const * const, bigint const * const);
status bigint_minus(bigint * const, bigint const * const, bigint const * const);
status bigint_mul(bigint * const, bigint const * const, bigint const * const);
status bigint_div(bigint * const, bigint const * const, bigint const * const);

// common functions
void print_err()
{
	printf("[error]");
}

// parser functions
lex create_lex()
{
	lex res;
	res.size = 0;
	res.type = EMPTY_LEX;
	res.val = NULL;
	return res;
}

void free_lex(lex * const l)
{
	if (l->val) {
		free(l->val);
		l->val = NULL;
		l->size = 0;
		l->type = EMPTY_LEX;
	}
	return;
}

status read_lex(lex * const l)
{
	if (!l)
		return RT_ERROR;
	free_lex(l);
	*l = create_lex();
	char cur;
	while (true)
	{
		cur = getchar();
		if (iswdigit(cur))
		{
			l->type = NUMBER;
			ungetc(cur, stdin);
			while (iswdigit(cur = getchar())) {
				l->val = (char *)realloc(l->val, ++l->size * sizeof(char));
				if (!l->val) {
					free_lex(l);
					return MEMORY_ERROR;
				}
				l->val[l->size - 1] = cur;
			}
			//l->val = (char *)realloc(l->val, ++l->size * sizeof(char));
			//if (!l->val) {
				//free_lex(l);
				//return MEMORY_ERROR;
			//}
			//l->val[l->size - 1] = '\0';
			ungetc(cur, stdin);
			return SUCCESS;
		}
		else if (cur == '(') {
			l->type = OPEN;
			return SUCCESS;
		}
		else if (cur == ')') {
			l->type = CLOSE;
			return SUCCESS;
		}
		else if (cur == '+') {
			l->type = PLUS;
			return SUCCESS;
		}
		else if (cur == '-') {
			l->type = MINUS;
			return SUCCESS;
		}
		else if (cur == '*') {
			l->type = MUL;
			return SUCCESS;
		}
		else if (cur == '/') {
			l->type = DIV;
			return SUCCESS;
		}
		else if (cur == EOF || cur == '\n') {
			l->type = EOF_LEX;
			return SUCCESS;
		}
		else if (!iswspace(cur) /*&& cur != '\n'*/) { // TODO: check if \n is a bad symbol (while tests)
			free_lex(l);
			return LEX_ERROR;
		}
	}

	return SUCCESS;
}

void print_lex(lex const * const l)
{
	switch (l->type) {
	case EMPTY_LEX:
		printf("EMPTY_LEX");
		break;
	case EOF_LEX:
		printf("EOF_LEX");
		break;
	case OPEN:
		printf("(");
		break;
	case CLOSE:
		printf(")");
		break;
	case PLUS:
		printf("+");
		break;
	case MINUS:
		printf("-");
		break;
	case MUL:
		printf("*");
		break;
	case DIV:
		printf("/");
		break;
	case NUMBER:
		for (size_t i = 0; i < l->size; ++i)
			printf("%c", l->val[i]);
		break;
	default:
		printf("Unknown lexeme!");
	}
	return;
}

status calc_expression(lex * const lexeme, bigint * const res)
{
	status st;
	free_bigint(res);
	free_lex(lexeme);
	// read first lexeme
	if ((st = read_lex(lexeme)) != SUCCESS)
		return st;
	if (lexeme->type == EOF_LEX)
		return SUCCESS;
	if ((st = expression(lexeme, res)) != SUCCESS) {
		free_lex(lexeme);
		free_bigint(res);
		return st;
	}
	if (lexeme->type != EOF_LEX) {
		free_lex(lexeme);
		free_bigint(res);
		return SYNT_ERROR;
	}
	free_lex(lexeme);
	return SUCCESS;
}

status expression(lex * const lexeme, bigint * const res)
{
	status st;
	if ((st = term(lexeme, res)) != SUCCESS)
		return st;
	while (lexeme->type == PLUS || lexeme->type == MINUS) {
		lex_type op_type = lexeme->type;
		if ((st = read_lex(lexeme)) != SUCCESS)
			return st;
		bigint tmp = create_bigint();
		if ((st = bigint_copy(&tmp, res)) != SUCCESS)
			return st;
		if ((st = term(lexeme, res)) != SUCCESS) {
			free_bigint(&tmp);
			return st;
		}
		bigint tmp_op_res = create_bigint();
		if (op_type == PLUS) {
			if ((st = bigint_plus(&tmp_op_res, &tmp, res)) != SUCCESS) {
				free_bigint(&tmp);
				return st;
			}
		}
		else {
			if ((st = bigint_minus(&tmp_op_res, &tmp, res)) != SUCCESS) {
				free_bigint(&tmp);
				return st;
			}
		}
		free_bigint(&tmp);
		if ((st = bigint_copy(res, &tmp_op_res)) != SUCCESS) {
			free_bigint(&tmp_op_res);
			return st;
		}
		free_bigint(&tmp_op_res);
	}
	return SUCCESS;
}

status term(lex * const lexeme, bigint * const res)
{
	status st;
	if ((st = factor(lexeme, res)) != SUCCESS)
		return st;
	while (lexeme->type == MUL || lexeme->type == DIV) {
		lex_type op_type = lexeme->type;
		if ((st = read_lex(lexeme)) != SUCCESS)
			return st;
		bigint tmp = create_bigint();
		if ((st = bigint_copy(&tmp, res)) != SUCCESS)
			return st;
		if ((st = factor(lexeme, res)) != SUCCESS) {
			free_bigint(&tmp);
			return st;
		}
		bigint tmp_op_res = create_bigint();
		if (op_type == MUL) {
			if ((st = bigint_mul(&tmp_op_res, &tmp, res)) != SUCCESS) {
				free_bigint(&tmp);
				return st;
			}
		}
		else {
			if ((st = bigint_div(&tmp_op_res, &tmp, res)) != SUCCESS) {
				free_bigint(&tmp);
				return st;
			}
		}
		free_bigint(&tmp);
		if ((st = bigint_copy(res, &tmp_op_res)) != SUCCESS) {
			free_bigint(&tmp_op_res);
			return st;
		}
		free_bigint(&tmp_op_res);
	}
	return SUCCESS;
}

status factor(lex * const lexeme, bigint * const res)
{
	status st;
	if (lexeme->type == MINUS)
	{
		if ((st = read_lex(lexeme)) != SUCCESS)
			return st;
		if ((st = factor(lexeme, res)) != SUCCESS)
			return st;
		res->is_negative = !res->is_negative;
	}
	else if (lexeme->type == OPEN) {
		if ((st = read_lex(lexeme)) != SUCCESS)
			return st;
		if ((st = expression(lexeme, res)) != SUCCESS)
			return st;
		if (lexeme->type != CLOSE)
			return SYNT_ERROR;
		// read next lexeme
		if ((st = read_lex(lexeme)) != SUCCESS)
			return st;
	}
	else if (lexeme->type == NUMBER) {
		if ((st = bigint_from_string(res, lexeme->val, lexeme->size)) != SUCCESS)
			return st;
		// read next lexeme
		if ((st = read_lex(lexeme)) != SUCCESS)
			return st;
	}
	else
		return SYNT_ERROR;
	return SUCCESS;
}


// long arithmetic functions
bigint create_bigint()
{
	bigint res;
	res.digits_size = 0;
	res.is_negative = false;
	res.digits = NULL;
	return res;
}

void free_bigint(bigint * const n)
{
	if (n->digits) {
		free(n->digits);
		n->digits = NULL;
		n->digits_size = 0;
		n->is_negative = false;
	}
	return;
}

void print_bigint(bigint const * const n)
{
	if (is_bigint_zero(n)) {
		printf("0");
		return;
	}
	if (n->is_negative)
		printf("-");
	size_t i = n->digits_size;
	while (i-- > 0)
		if (i == n->digits_size - 1)
			printf("%" PRIdMAX "", (intmax_t)n->digits[i]);
		else
			printf("%04" PRIdMAX "", (intmax_t)n->digits[i]);
	return;
}

bool is_bigint_zero(bigint const * const n)
{
	return n->digits_size == 1 && !n->digits[0];
}

bool is_left_bigint_ge_abs(bigint const * const left, bigint const * const right)
{
	if (left->digits_size > right->digits_size)
		return true;
	if (left->digits_size < right->digits_size)
		return false;
	size_t i = left->digits_size;
	for (; i > 0; --i) {
		if (left->digits[i - 1] > right->digits[i - 1])
			return true;
		if (left->digits[i - 1] < right->digits[i - 1])
			return false;
	}
	return true;
}

status n_digits_from_string(bigint_digit * const dst, char const * const str, size_t const n)
{
	*dst = 0;
	for (size_t i = 0; i < n; ++i) {
		*dst *= 10;
		*dst += str[i] - '0';
	}
	return SUCCESS;
}

status bigint_from_string(bigint * const res, char const * const str, size_t const size)
{
	free_bigint(res);
	// delete leading zeros
	intptr_t start = 0;
	while ((size_t)start < size && str[start++] == '0');
	start--;
	intptr_t no_lead_zeros_size = size - start;
	// allocate memory for bigint digits array
	size_t digits_size = no_lead_zeros_size / BIGINT_BASE_DEC;
	digits_size += (no_lead_zeros_size%BIGINT_BASE_DEC) ? 1 : 0;
	res->digits = (bigint_digit *)malloc(digits_size * sizeof(bigint_digit));
	if (!res->digits)
		return MEMORY_ERROR;
	res->digits_size = digits_size;
	// fill digits with values
	digits_size = 0;
	for (intptr_t i = size - 1; i >= start; i -= BIGINT_BASE_DEC) {
		int8_t n_digits = (i - BIGINT_BASE_DEC < start) ? i + 1 - start : BIGINT_BASE_DEC;
		intptr_t start_idx = (i - BIGINT_BASE_DEC < start) ? start : i + 1 - BIGINT_BASE_DEC;
		status st = n_digits_from_string(&(res->digits[digits_size++]), str + start_idx, (size_t)n_digits);
		if (st != SUCCESS) {
			free_bigint(res);
			return st;
		}
	}
	return SUCCESS;
}

status bigint_copy(bigint * const dst, bigint const * const src)
{
	free_bigint(dst);
	dst->digits = (bigint_digit *)malloc(src->digits_size * sizeof(bigint_digit));
	if (!dst->digits)
		return MEMORY_ERROR;
	dst->digits_size = src->digits_size;
	dst->is_negative = src->is_negative;
	memcpy(dst->digits, src->digits, src->digits_size * sizeof(bigint_digit));
	return SUCCESS;
}

status bigint_abs(bigint * const dst, bigint const * const src)
{
	status st;
	free_bigint(dst);
	if ((st = bigint_copy(dst, src)) != SUCCESS)
		return st;
	dst->is_negative = false;
	return SUCCESS;
}

status bigint_delete_lead_zeros(bigint * const res)
{
	while (res->digits_size > 1 && !res->digits[res->digits_size - 1])
		res->digits_size--;
	return SUCCESS;
}

status bigint_insert(bigint * const res, size_t const pos, bigint_digit const digit)
{
	status st;
	if (pos > res->digits_size)
		return RT_ERROR;
	bigint_digit * tmp_digits = (bigint_digit *)realloc(res->digits, (res->digits_size + 1) * sizeof(bigint_digit));
	if (!tmp_digits)
		return MEMORY_ERROR;
	res->digits = tmp_digits;
	res->digits_size++;

	tmp_digits = (bigint_digit *)memmove(res->digits + pos + 1, res->digits + pos,
		(res->digits_size - pos - 1) * sizeof(bigint_digit));
	if (!tmp_digits)
		return MEMORY_ERROR;
	res->digits[pos] = digit;
	if ((st = bigint_delete_lead_zeros(res)) != SUCCESS)
		return st;
	return SUCCESS;
}

status bigint_push_back(bigint * const res, bigint_digit const digit)
{
	status st;
	bigint_digit * tmp_digits = (bigint_digit *)realloc(res->digits, (res->digits_size + 1) * sizeof(bigint_digit));
	if (!tmp_digits)
		return MEMORY_ERROR;
	res->digits = tmp_digits;
	res->digits[res->digits_size++] = digit;
	if ((st = bigint_delete_lead_zeros(res)) != SUCCESS)
		return st;
	return SUCCESS;
}

status bigint_push_front(bigint * const res, bigint_digit const digit)
{
	status st;
	if ((st = bigint_insert(res, 0, digit)) != SUCCESS)
		return st;
	return SUCCESS;
}

status bigint_nonzero_nonneg_plus(bigint * const res, bigint const * const bigger, bigint const * const smaller)
{
	// allocate memory for res and choose bigger and smaller terms
	status st;
	free_bigint(res);
	if ((st = bigint_copy(res, bigger)) != SUCCESS)
		return st;
	// add digits
	bool carry = 0;
	for (size_t i = 0; i < bigger->digits_size || carry; ++i) {
		if (i == res->digits_size) {
			if ((st = bigint_push_back(res, 0)) != SUCCESS) {
				free_bigint(res);
				return st;
			}
		}
		res->digits[i] += carry + (i < smaller->digits_size ? smaller->digits[i] : 0);
		if ((carry = (res->digits[i] >= BIGINT_BASE)))
			res->digits[i] -= BIGINT_BASE;
	}
	return SUCCESS;
}

status bigint_nonzero_nonneg_minus(bigint * const res, bigint const * const bigger, bigint const * const smaller)
{
	// allocate memory for res and choose bigger and smaller terms
	status st;
	free_bigint(res);
	if ((st = bigint_copy(res, bigger)) != SUCCESS)
		return st;
	// subtract digits
	bool carry = 0;
	for (size_t i = 0; i < smaller->digits_size || carry; ++i) {
		res->digits[i] -= carry + (i < smaller->digits_size ? smaller->digits[i] : 0);
		carry = res->digits[i] < 0;
		if (carry)
			res->digits[i] += BIGINT_BASE;
	}
	// delete leading zeros
	if ((st = bigint_delete_lead_zeros(res)) != SUCCESS) {
		free_bigint(res);
		return st;
	}

	return SUCCESS;
}

status bigint_nonzero_nonneg_mul(bigint * const res, bigint const * const bigger, bigint const * const smaller)
{
	// allocate memory for res and fill it with zeros
	status st;
	free_bigint(res);
	res->digits = (bigint_digit *)malloc((bigger->digits_size + smaller->digits_size) * sizeof(bigint_digit));
	if (!res->digits)
		return MEMORY_ERROR;
	res->digits_size = bigger->digits_size + smaller->digits_size;
	memset(res->digits, 0, res->digits_size * sizeof(bigint_digit));

	// mul digits
	bigint_digit carry = 0;
	bigint_digit tmp;
	for (size_t i = 0; i < smaller->digits_size; ++i) {
		carry = 0;
		for (size_t j = 0; j < bigger->digits_size || carry; ++j) {
			bigint_digit bigger_digits = (j < bigger->digits_size ? bigger->digits[j] : 0);
			tmp = res->digits[i + j] + smaller->digits[i] * bigger_digits + carry;
			res->digits[i + j] = tmp % BIGINT_BASE;
			carry = tmp / BIGINT_BASE;
		}
	}
	// delete leading zeros
	if ((st = bigint_delete_lead_zeros(res)) != SUCCESS) {
		free_bigint(res);
		return st;
	}

	return SUCCESS;
}

status bigint_nonzero_nonneg_div(bigint * const res, bigint const * const bigger, bigint const * const smaller)
{
	status st;
	free_bigint(res);
	// allocate memory for res and fill it with zeros
	res->digits = (bigint_digit *)malloc((bigger->digits_size) * sizeof(bigint_digit));
	if (!res->digits)
		return MEMORY_ERROR;
	res->digits_size = bigger->digits_size;
	memset(res->digits, 0, res->digits_size * sizeof(bigint_digit));

	// div digits
	bigint cur_bigint = create_bigint();
	bigint tmp_bigint = create_bigint();
	bigint tmp2_bigint = create_bigint();
	bigint m_bigint = create_bigint();
	bigint x_bigint = create_bigint();
	for (size_t cur_digit_num = bigger->digits_size; cur_digit_num > 0; --cur_digit_num) {
		//  update cur_bigint with new digits
		if ((st = bigint_push_front(&cur_bigint, bigger->digits[cur_digit_num - 1])) != SUCCESS) {
			free_bigint(&cur_bigint);
			return st;
		}
		// if cur_bigint >= smaller, use binary search to find max x: smaller * x <= cur_bigint
		if (is_left_bigint_ge_abs(&cur_bigint, smaller)) {
			bigint_digit x = 0;
			bigint_digit left = 0;
			bigint_digit right = BIGINT_BASE;
			while (left <= right) {
				bigint_digit middle = (left + right) >> 1;
				if ((st = bigint_push_back(&m_bigint, middle)) != SUCCESS) {
					free_bigint(res);
					free_bigint(&cur_bigint);
					return st;
				}
				if ((st = bigint_mul(&tmp_bigint, smaller, &m_bigint)) != SUCCESS) {
					free_bigint(res);
					free_bigint(&cur_bigint);
					free_bigint(&m_bigint);
					return st;
				}
				if (is_left_bigint_ge_abs(&cur_bigint, &tmp_bigint)) {
					x = middle;
					left = middle + 1;
				}
				else
					right = middle - 1;
				free_bigint(&tmp_bigint);
				free_bigint(&m_bigint);
			}
			// store x to answer
			res->digits[cur_digit_num - 1] = x;
			// update cur_bigint
			if ((st = bigint_push_back(&x_bigint, x)) != SUCCESS) {
				free_bigint(res);
				free_bigint(&cur_bigint);
				return st;
			}
			if ((st = bigint_mul(&tmp_bigint, smaller, &x_bigint)) != SUCCESS) {
				free_bigint(res);
				free_bigint(&cur_bigint);
				free_bigint(&x_bigint);
				return st;
			}
			free_bigint(&x_bigint);
			if ((st = bigint_minus(&tmp2_bigint, &cur_bigint, &tmp_bigint)) != SUCCESS) {
				free_bigint(res);
				free_bigint(&cur_bigint);
				free_bigint(&tmp_bigint);
				return st;
			}
			free_bigint(&tmp_bigint);
			if ((st = bigint_copy(&cur_bigint, &tmp2_bigint)) != SUCCESS) {
				free_bigint(res);
				free_bigint(&cur_bigint);
				free_bigint(&tmp2_bigint);
				return st;
			}
			free_bigint(&tmp2_bigint);
		}
	}
	free_bigint(&cur_bigint);
	// delete leading zeros
	if ((st = bigint_delete_lead_zeros(res)) != SUCCESS) {
		free_bigint(res);
		return st;
	}
	return SUCCESS;
}

status bigint_plus(bigint * const res, bigint const * const left, bigint const * const right)
{
	status st;
	// if one of terms is zero
	if (is_bigint_zero(left)) {
		if ((st = bigint_copy(res, right)) != SUCCESS)
			return st;
		return SUCCESS;
	}
	if (is_bigint_zero(right)) {
		if ((st = bigint_copy(res, left)) != SUCCESS)
			return st;
		return SUCCESS;
	}
	// for non-zero terms
	// determine bigger and smaller terms
	bool left_is_ge = is_left_bigint_ge_abs(left, right);
	bigint const * const bigger = left_is_ge ? left : right;
	bigint const * const smaller = left_is_ge ? right : left;

	if (bigger->is_negative != smaller->is_negative) {
		if ((st = bigint_nonzero_nonneg_minus(res, bigger, smaller)) != SUCCESS)
			return st;
	}
	else {
		if ((st = bigint_nonzero_nonneg_plus(res, bigger, smaller)) != SUCCESS)
			return st;
	}
	res->is_negative = bigger->is_negative;

	return SUCCESS;
}

status bigint_minus(bigint * const res, bigint const * const left, bigint const * const right)
{
	status st;
	// if one of terms is zero
	if (is_bigint_zero(left)) {
		if ((st = bigint_copy(res, right)) != SUCCESS)
			return st;
		res->is_negative = !right->is_negative;
		return SUCCESS;
	}
	if (is_bigint_zero(right)) {
		if ((st = bigint_copy(res, left)) != SUCCESS)
			return st;
		return SUCCESS;
	}
	// for non-zero terms
	// determine bigger and smaller terms
	bool left_is_ge = is_left_bigint_ge_abs(left, right);
	bigint const * const bigger = left_is_ge ? left : right;
	bigint const * const smaller = left_is_ge ? right : left;
	if (bigger->is_negative != smaller->is_negative) {
		if ((st = bigint_nonzero_nonneg_plus(res, bigger, smaller)) != SUCCESS)
			return st;
	}
	else {
		if ((st = bigint_nonzero_nonneg_minus(res, bigger, smaller)) != SUCCESS)
			return st;
	}
	res->is_negative = bigger->is_negative == left_is_ge;

	return SUCCESS;
}

status bigint_mul(bigint * const res, bigint const * const left, bigint const * const right)
{
	status st;
	// if one of factors is zero
	if (is_bigint_zero(left)) {
		if ((st = bigint_copy(res, left)) != SUCCESS)
			return st;
		return SUCCESS;
	}
	if (is_bigint_zero(right)) {
		if ((st = bigint_copy(res, right)) != SUCCESS)
			return st;
		return SUCCESS;
	}
	// for non-zero factors
	// determine bigger and smaller terms
	bool left_is_ge = is_left_bigint_ge_abs(left, right);
	bigint const * const bigger = left_is_ge ? left : right;
	bigint const * const smaller = left_is_ge ? right : left;
	if ((st = bigint_nonzero_nonneg_mul(res, bigger, smaller)) != SUCCESS)
		return st;
	res->is_negative = bigger->is_negative != smaller->is_negative;

	return SUCCESS;
}

status bigint_div(bigint * const res, bigint const * const left, bigint const * const right)
{
	status st;
	// if one of factors is zero
	if (is_bigint_zero(left) && !is_bigint_zero(right)) {
		if ((st = bigint_copy(res, left)) != SUCCESS)
			return st;
		return SUCCESS;
	}
	if (is_bigint_zero(right)) {
		return RT_ERROR;
	}
	// for non-zero factors
	if (!is_left_bigint_ge_abs(left, right)) {
		free_bigint(res);
		if ((st = bigint_push_back(res, 0)) != SUCCESS)
			return st;
	}
	else {
		// work with abs of factors
		bigint left_abs = create_bigint();
		bigint right_abs = create_bigint();
		if ((st = bigint_abs(&left_abs, left)) != SUCCESS)
			return st;
		if ((st = bigint_abs(&right_abs, right)) != SUCCESS) {
			free_bigint(&left_abs);
			return st;
		}
		if ((st = bigint_nonzero_nonneg_div(res, &left_abs, &right_abs)) != SUCCESS)
			return st;
		res->is_negative = left->is_negative != right->is_negative;
		free_bigint(&left_abs);
		free_bigint(&right_abs);
	}

	return SUCCESS;
}


int main(int argc, char ** argv)
{
	status st;
	lex lexeme = create_lex();
	bigint res = create_bigint();
	if ((st = calc_expression(&lexeme, &res) != SUCCESS) ){
		print_err();
		return SUCCESS;
	}
	print_bigint(&res);
	free_bigint(&res);
	free_lex(&lexeme);
	return SUCCESS;
}