#ifndef BIO_UNARY_COMPOSE_H_
#define BIO_UNARY_COMPOSE_H_

#include "bio/defs.h"

#include <functional>

BIO_NS_START




template< typename Op1, typename Op2 >
class unary_compose
{
public:
	typedef typename Op1::result_type result_type;

protected:
	Op1 op1;
	Op2 op2;

public:
	unary_compose(
		const Op1 & op1 = Op1(),
		const Op2 & op2 = Op2() ) 
		: op1( op1 )
		, op2( op2 )
	{
	}

	template< typename argument_type >
	result_type &
	operator()( const argument_type& x ) const
	{
		return op1( op2( x ) );
	}

	template< typename argument_type >
	result_type &
	operator()( const argument_type& x )
	{
		return op1( op2( x ) );
	}
};


template< typename Op1, typename Op2 >
unary_compose< Op1, Op2 >
compose(
	const Op1 & op1,
	const Op2 & op2 )
{
	return unary_compose< Op1, Op2 >( op1, op2 );

}


BIO_NS_END

#endif //BIO_UNARY_COMPOSE_H_
