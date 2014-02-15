
template<typename Type>
class myMat {

public:
	//constructor and destructor
     myMat();
     myMat( const myMat<Type> &A );
     myMat( int rows, int columns, int slices=0 ); //const Type &x = Type(0)
    // myMat( int rows, int columns,const Type *v);
     ~myMat();

        // assignments
    myMat<Type>& operator=( const myMat<Type> &A );
    myMat<Type>& operator=( const Type &x );

        // accessors
    Type operator[]( int i );
    const Type operator[]( int i ) const;
    Type& operator()( int row, int column,int slice=0 );
    const Type& operator()( int row, int column, int slice=0) const;

        // type conversion
    operator Type*();
    operator const Type*() const;

        // others
    long size() const;
    int dim( int dimension ) const;
    int rows() const;
    int cols() const;
    int slices() const;

private:
	int dimension;
	int nRow;
	int nColumn;
	int nSlice;
	int nTotal;
	Type *data;


}