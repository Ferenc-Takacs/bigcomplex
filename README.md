Use python command:

from complex import *

The precision set or get by __getcontext().prec__ variable.

Content:
```bash
def binom(x,y): Binomial coeficient.
					y can int type,
					x can int, float, decimal, complex, Real, Complex types

def pi() :
    	return Complex(Real.pi())

def e():
    	return Complex(Real.e())

def j():
    	return Complex(1j)

class Real(Decimal) : #extensions for decimal
	  def pi() :
	  def e():
      def cos(self):
      def sin(self):
      def cos(self):
      def tan(self):
      def atan2(y,x):
      def atan(self):
      def asin(self):
      def acos(self):
      def actg(self):
      def grad(self):
      def rad(self):
      def sqrt(self):

class Complex(object):
	    def __init__(self,real,imag=0):
    	def pi() :
    	def e():
    	def j():
    	def conj(self):
    	def re(self):
    	def im(self):
    	def floor(self):
    	def ceil(self):
    	def round(self):
    	def exp(self):
    	def phase(self):
    	def polar(self):
	    def rect(self,r,phi):
     	def cos(self):
     	def sin(self):
     	def tan(self):
     	def acos(self):
     	def asin(self):
     	def atan(self):
     	def sinh(self):
     	def cosh(self):
     	def tanh(self):
     	def acosh(self):
     	def asinh(self):
     	def atanh(self):
     	def ln(self):
     	def __pow__(self,other):    # Complex ** other
	    def __rpow__(self,other):	# other ** Complex
     	def __str__(self):			# to string
     	def __repr__(self):			# debug string
     	def __abs__(self):			# abs( Complex )
     	def abs(self):				# Complex.abs()
     	def __add__(self,other):	# Complex + other
     	def __radd__(self,other):	# other + Complex
     	def __sub__(self,other):	# Complex - other
     	def __rsub__(self,other):	# other - Complex
     	def __mul__(self,other):	# Complex * other
     	def __rmul__(self,other):	# other * Complex
     	def __floordiv__(self,other):	# Complex.real() // other
     	def __rfloordiv__(self,other):	# Complex(other) // Complex
     	def __truediv__(self,other):	# Complex / other
     	def __rtruerdiv__(self,other):	# other / Complex
     	def __neg__(self):			# -Complex
     	def __pos__(self):			# normalizale imagine and real part together
     	def __inverse__(self):		# 1 / Complex
     	def __mod__(self,other):	# Complex % other   #  a + b * ciel(-a/b)
     	def __complex__(self):		# complex(Complex)
     	def __int__(self):			# int(Complex)
     	def __float__(self):		# float(Complex)
     	def __lt__(self,other):		# Complex < other
     	def __le__(self,other):		# Complex <= other
     	def __eq__(self,other):		# Complex == other
     	def __ne__(self,other):		# Complex != other
     	def __gt__(self,other):		# Complex > other
     	def __ge__(self,other):		# Complex >= other
     	def is_finite(self):
     	def is_infinite(self):
     	def is_nan(self):
     	def is_signed(self):
     	def is_zero(self):
     	def sqrt(self):
     	def to_eng_string(self):
```
