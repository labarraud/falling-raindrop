
#ifndef LINALG_FILE_VECTOR_CXX

#include "Vector.hxx"

namespace linalg
{
  
  //! Default constructor : empty vector
  template<class T, class Allocator>
  inline Vector<T, Allocator>::Vector()
  {
    this->m_ = 0;
    this->data_ = NULL;
  }
  
  
  //! Constructor with the size of the vector
  template<class T, class Allocator>
  inline Vector<T, Allocator>::Vector(size_t n)
  {
#ifdef LINALG_DEBUG
    try
      {
#endif
	this->m_ = n;
	this->data_ = Allocator::allocate(n);
	
#ifdef LINALG_DEBUG
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    
    if (this->data_ == NULL)
      this->m_ = 0;
    
    if (this->data_ == NULL && n != 0)
      throw NoMemory("Vector::Vector(size_t)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_string(n*sizeof(T)) + " bytes ("
		     + to_string(n) + " elements).");
#endif
  }
  
  
  //! Copy constructor
  template<class T, class Allocator>
  inline Vector<T, Allocator>::Vector(const Vector<T, Allocator>& V)
  {
#ifdef LINALG_DEBUG
    try
      {
#endif
	this->m_ = V.GetSize();
	this->data_ = Allocator::allocate(V.GetSize());

#ifdef LINALG_DEBUG
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    
    if (this->data_ == NULL)
      this->m_ = 0;
    
    if (this->data_ == NULL && V.GetSize() != 0)
      throw NoMemory("Vector::Vector(Vector&)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_string(V.GetSize()*sizeof(T)) + " bytes ("
		     + to_string(V.GetSize()) + " elements).");
#endif
    
    Allocator::memorycpy(this->data_, V.GetData(), V.GetSize());
  }
  
  
  //! Destructor
  template<class T, class Allocator>
  inline Vector<T, Allocator>::~Vector()
  {
    Clear();
  }
  
  
  //! Clears the vector
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::Clear()
  {
#ifdef LINALG_DEBUG
    try
      {
#endif
	
	if (data_ != NULL)
	  {
	    Allocator::deallocate(data_, m_);
	    m_ = 0;
	    data_ = NULL;
	  }
	
#ifdef LINALG_DEBUG
      }
    catch (...)
      {
	m_ = 0;
	data_ = NULL;
      }
#endif
  }
  

  //! Returns the size of the vector (as an int)
  template<class T, class Allocator>
  inline int Vector<T, Allocator>::GetM() const
  {
    return this->m_;
  }
  
  
  //! Returns the size of the vector
  template<class T, class Allocator>
  inline size_t Vector<T, Allocator>::GetSize() const
  {
    return this->m_;
  }
  
  
  //! Returns the C pointer associated with the vector
  template<class T, class Allocator>
  inline T* Vector<T, Allocator>::GetData() const
  {
    return this->data_;
  }
  
  
  //! Returns the C pointer associated with the vector
  template<class T, class Allocator>
  inline void* Vector<T, Allocator>::GetDataVoid() const
  {
    return reinterpret_cast<void*>(this->data_);
  }
  
  
  //! Changes the size of the vector (previous elements are not kept)
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::Reallocate(size_t n)
  {
#ifdef LINALG_DEBUG
    try
      {
#endif
	
	if (n < this->m_)
	  {
	    Clear();
	    
	    this->m_ = n;
	    this->data_ = Allocator::allocate(n);
	  }      
	else if (n > this->m_)
	  {
	    this->m_ = n;
	    
	    this->data_ =
	      reinterpret_cast<T*>(Allocator::
				   reallocate(this->data_, n));
	    
	  }
	
#ifdef LINALG_DEBUG
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    
    if (this->data_ == NULL)
      this->m_ = 0;
    
    if (this->data_ == NULL && n != 0)
      throw NoMemory("Vector::Reallocate(size_t)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_string(n*sizeof(T)) + " bytes ("
		     + to_string(n) + " elements).");
#endif
  }
  

  //! Changes the size of the vector (previous elements are kept)  
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::Resize(size_t n)
  {
#ifdef LINALG_DEBUG
    try
      {
#endif
	
	this->data_ =
	  reinterpret_cast<T*>(Allocator::
			       resize(this->data_, this->m_, n));
	
	this->m_ = n;
#ifdef LINALG_DEBUG
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    
    if (this->data_ == NULL)
      this->m_ = 0;
    
    if (this->data_ == NULL && n != 0)
      throw NoMemory("Vector::Resize(size_t)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_string(n*sizeof(T)) + " bytes ("
		     + to_string(n) + " elements).");
    
#endif
  }
  
  
  //! Sets the size of the vector and the associated pointer
  /*!
    This function is a low-level function and should be used cautiously
   */
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::SetData(size_t n, T* data)
  {
    Clear();
    data_ = data;
    m_ = n;
  }
  

  //! Sets the C pointer to 0
  /*!
    This function is a low-level function and should be used cautiously
   */  
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::Nullify()
  {
    m_ = 0;
    data_ = NULL;
  }
  
  
  //! Returns access to the element i
  template<class T, class Allocator>
  inline T& Vector<T, Allocator>::operator()(size_t i)
  {
#ifdef LINALG_DEBUG
    if (i >= this->m_)
      throw WrongIndex("Vector::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_string(int(this->m_)-1) + "], but is equal to "
		       + to_string(i) + ".");
#endif
    
    return this->data_[i];
  }
  
  
  //! Returns access to the element i
  template<class T, class Allocator>
  inline const T& Vector<T, Allocator>::operator()(size_t i) const
  {
#ifdef LINALG_DEBUG
    if (i >= this->m_)
      throw WrongIndex("Vector::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_string(int(this->m_)-1) + "], but is equal to "
		       + to_string(i) + ".");
#endif
    
    return this->data_[i];
  }
    
  
  //! Operator =, *this = X
  template<class T, class Allocator>
  inline Vector<T, Allocator>& Vector<T, Allocator>::operator=(const Vector<T, Allocator>& X)
  {
    this->Reallocate(X.GetSize());
    
    Allocator::memorycpy(this->data_, X.GetData(), this->m_);
    return *this;
  }
  

  //! Multiplies the vector by a scalar
  template<class T, class Allocator>
  Vector<T, Allocator>& Vector<T, Allocator>::operator*=(const T& alpha)
  {
    for (size_t i = 0; i < this->GetSize(); i++)
      this->data_[i] *= alpha;
    
    return *this;
  }
  
  
  //! Appends x at the end of the vector
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::PushBack(const T& x)
  {
    Resize(this->m_+1);
    this->data_[this->m_-1] = x;
  }
  
  
  //! Appends x at the end of the vector
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::PushBack(const Vector<T, Allocator>& X)
  {
    int Nold = this->m_;
    Resize(this->m_ + X.GetSize());
    for (size_t i = 0; i < X.GetSize(); i++)
      this->data_[Nold+i] = X(i);
  }
  
  
  //! Sets the vector with 0 (valid for basic types)
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::Zero()
  {
    Allocator::memoryset(this->data_, char(0),
			 this->GetSize() * sizeof(T));
  }
  

  //! Fills the vector with the same value x
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::Fill(const T& x)
  {
    for (size_t i = 0; i < this->m_; i++)
      data_[i] = x;
  }
  

  //! Fills the vector randomly
  template<class T, class Allocator>
  inline void Vector<T, Allocator>::FillRand()
  {
    for (size_t i = 0; i < this->m_; i++)
      GetRand(data_[i]);
  }
  
  
  //! Writes the vector in a file.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector::Write(string FileName, "
                    "bool with_size)",
		    string("Unable to open file \"") + FileName + "\".");

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector::Write(ostream& FileStream, "
                    "bool with_size)",
                    "The stream is not ready.");

    if (with_size)
      FileStream.write(reinterpret_cast<char*>(const_cast<size_t*>(&this->m_)),
		       sizeof(size_t));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * sizeof(T));

    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector::Write(ostream& FileStream, "
                    "bool with_size)",
                    "Output operation failed.");

  }
  
  
  //! Writes the content of the vector in a text file
  template<class T, class Allocator>
  void Vector<T, Allocator>::WriteText(const string& FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
    
    this->WriteText(FileStream);
    
    // end of line to finish the file
    FileStream << '\n';
    
    FileStream.close();    
  }
  
  
  //! Writes the content of the vector in an output stream
  template<class T, class Allocator>
  void Vector<T, Allocator>::WriteText(ostream& FileStream) const
  {
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
    
    for (size_t i = 0; i < this->GetSize(); i++)
      FileStream << this->data_[i] << '\t';
    
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector::WriteText(ostream& FileStream)",
                    "Output operation failed.");
  }
    
  
  //! Reads the contents of the vector in a text file
  template<class T, class Allocator>
  void Vector<T, Allocator>::ReadText(const string& FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");

    this->ReadText(FileStream);

    FileStream.close();

  }
  

  //! Reads the contents of the vector in an input stream
  template<class T, class Allocator>
  void Vector<T, Allocator>::ReadText(istream& FileStream)
  {
    // Previous values of the vector are cleared.
    Clear();

    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::ReadText(istream& FileStream)",
                    "The stream is not ready.");

    T entry;
    size_t number_element = 0;
    while (!FileStream.eof())
      {
	// Reads a new entry.
	FileStream >> entry;
	if (FileStream.fail())
	  break;
	else
	  {
	    number_element++;

	    // If needed, resizes the vector. Its size is already doubled so
	    // that the vector should be resized a limited number of times.
	    if (number_element > this->m_)
	      this->Resize(2 * number_element);

	    this->data_[number_element - 1] = entry;
	  }
      }
    
    // Resizes to the actual size.
    if (number_element > 0)
      this->Resize(number_element);
    else
      this->Clear();

  }


  template<class T>
  inline void GetRand(T& x)
  {
    x = T(rand()) / RAND_MAX;
  }

  
  template<class T>
  void GetRand(complex<T>& x)
  {
    int p = rand()%3;
    if (p == 0)
      x = complex<T>(rand(), 0) / T(RAND_MAX);
    else if (p == 1)
      x = complex<T>(0, rand()) / T(RAND_MAX);
    else
      x = complex<T>(rand(), rand()) / T(RAND_MAX);
  }


  template<class T>
  T DotProd(const Vector<T>& x, Vector<T>& y)
  {
    T sum(0);
    for (int i = 0; i < x.GetM(); i++)
      sum += x(i)*y(i);
    
    return sum;
  }

  template<class T>
  T Norm2(const Vector<T>& x)
  {
    T sum(0);
    for (int i = 0; i < x.GetM(); i++)
      sum += x(i)*x(i);
    
    return sqrt(sum);
  }

  template<class T>
  T Norm2(const Vector<complex<T> >& x)
  {
    T sum(0);
    for (int i = 0; i < x.GetM(); i++)
      sum += real(x(i))*real(x(i)) + imag(x(i))*imag(x(i));
    
    return sqrt(sum);
  }
  
  template<class T>
  void Add(const T& alpha, const Vector<T>& x, Vector<T>& y)
  {
    for (int i = 0; i < y.GetM(); i++)
      y(i) += alpha*x(i);
  }

	void Add(double alpha, const Vector<Vector<double> >& x, Vector<Vector<double> >& y)
  {
    for (int i = 0; i < y.GetM(); i++)
      Add(alpha, x(i), y(i));
  }

  template<class T, class Allocator>
  ostream& operator<<(ostream& out, const Vector<T, Allocator>& V)
  {
    V.WriteText(out);
    return out;
  }
  
}

#define LINALG_FILE_VECTOR_CXX
#endif
