import math
from decimal import *
import threading

def binom(x,y):
	inttype = isinstance(x,int)
	integer = inttype or (isinstance(x,Real) and x.is_int())
	if not(isinstance(x,(int,float,complex,Real,Complex))) :
		print('binom() error. First arg must be int, float, complex, Real(Decimal) or Complex.')
		return 0
	if not(isinstance(y,int)) :
		print('binom() error. Second arg must be int.')
		return 0
	if y < 0 :
		return 0
	if integer and x >= 0 :
		if y > x:
			return 0
		if 2*y > x :
			y = x - y
	if y == 0 or y == x:
		return 1
	if y == 1 or y == x-1:
		return x
	oldprec = getcontext().prec
	ax = int(abs(x))
	getcontext().prec = max(2 + (y+ax)//3, oldprec)
	newprec = getcontext().prec
	if isinstance(x,(complex,Complex)) :
		if isinstance(x,complex) :
			x = Complex(x)
	else :
		x = Real(x)
	result = x
	x -= 1
	i = 2
	while i <= y:
		result = result * x
		result = result / i
		x -= 1
		i += 1
	getcontext().prec = oldprec
	if inttype :
		return int(result)
	return result

def pi() :
	return Complex(Real.pi())

def e():
	return Complex(Real.e())

def j():
	return Complex(1j)

_Real_pi = 3.1
_Real_pi_prec = 0
_Real_e = 2.7
_Real_e_prec = 0
_Real_add_prec = 2
_Real_lock = threading.Lock()

class add_prec() :
	def __init__(self) :
		pass
	def __enter__(self) :
		self.pi = _Real_pi
		self.pi_prec = _Real_pi_prec
		self.e = _Real_e
		self.e_prec = _Real_e_prec
		getcontext().prec += _Real_add_prec
		return self
	def __exit__(self, exc_type, exc_value, traceback) :
		getcontext().prec -= _Real_add_prec
		_Real_pi = self.pi
		_Real_pi_prec = self.pi_prec
		_Real_e = self.e
		_Real_e_prec = self.e_prec

"""
The high precision number class extension for Decimal class
"""
class Real(Decimal) :

	def pi() :
		with _Real_lock :
			with add_prec() as prec:
				if _Real_pi_prec != getcontext().prec :
					lasts, t, result, n, na, d, da = 0, Real(3), Real(3), 1, 0, 0, 24
					while result != lasts:
						lasts = result
						n, na = n+na, na+8
						d, da = d+da, da+32
						t = (t * n) / d
						result += t
					prec.pi = +result # unary plus applies the new precision
					prec.pi_prec = getcontext().prec
				return prec.pi

	def e():
		with _Real_lock :
			with add_prec() as prec:
				if _Real_e_prec != getcontext().prec :
					prec.e = +Decimal(1).exp()
					prec.e_prec = getcontext().prec
				return prec.e

	def cos(x):
		pi = Real.pi()
		with add_prec():
			x = getcontext().remainder(abs(x),pi*2)
			if x > pi :
				x = 2*pi - x
			sign1 = 1
			if x*2 > pi :
				sign1 = -1
				x = pi - x
			i, lasts, result, fact, num, sign = Real(0), Real(0), Real(1), Real(1), Real(1), Real(1)
			while result != lasts:
				lasts = result
				i += 2
				fact *= i * (i-1)
				num *= x * x
				sign *= -1
				result += num / fact * sign
		return Real.filt(result*sign1)

	def sin(x):
		sin = Real.cos(x - Real.pi()/2)
		return Real.filt(sin)

	def tan(x):
		return Real.sin(x) / Real.cos(x)

	def atan2(y,x) :
		with add_prec():
			x = Real(x)
			y = Real(y)
			if abs(x) >= abs(y) :
				atan = Real.atan(y / x)
				#Complex.testprint(atan,'atan2(y/x)')
				if x < 0 :
					if y < 0 :
						atan -= Real.pi()
					else :
						atan += Real.pi()
			else :
				atan = Real.pi()/2 - Real.atan(x / y)
				#Complex.testprint(atan,'atan2(x/y)')
				if y < 0 :
					atan -= Real.pi()
		return Real.filt(atan)

	def atan(x):
		with add_prec():
			x = Real(x)
			if x.is_infinite() :
				if x.as_tuple()[0] > 0 :
					return Real.pi()/2
				else :
					return Real.pi()/-2
			#Complex.testprint(x,'atan x')
			d = x / Real.sqrt( (x*x) + 1 )
		return Real.asin( d )

	def asin(x) :
		with add_prec():
			x = Real(x)
			sign = 1 if x >= 0 else -1
			if x == 0 : return Real(0)
			pi = Real.pi()
			if x == 1 : return pi*Real(0.5)
			if x == -1 : return pi*Real(-0.5)
			if x > 1 : return 1/0
			if x < -1 : return -1/0
			x = abs(x)
			result, last, tag, n = x, x+1, x, 1
			#Complex.testprint(result,'asin x')
			while last!=result :
				last = result
				tag *= n * x * x
				n += 2
				tag /= n * (n-1)
				result += tag
			# very rough result
			xe = Real.sin(result) # inverse interpolation
			err = x - xe
			limit = Real(10)**Real(-getcontext().prec)
			#print( limit, result, xe, ' limit, asin res x, sin')
			result2 = result + result*err
			while abs(err)>limit :
				if result2 > pi : result2 = pi
				if result2 < -pi : result2 = -pi
				#print(result2,err,' err')
				xe2 = Real.sin(result2)
				err = x-xe2
				if err == 0 : break
				if xe2-xe == 0 : break
				result3 = err*(result2-result)/(xe2-xe) + result2
				result = result2
				result2 = result3
				xe = xe2
		return Real.filt(result2*sign)

	def acos(x) :
		acos = Real.pi()*Real(0.5) - Real.asin(x)
		return Real.filt(acos)

	def actg(x) :
		return Real.pi()*Real(0.5) - Real.atan(x)

	#def ln(x): # exist in Decimal.decimal
	
	#def exp(x): # exist in Decimal.decimal
	
	def grad(x) :
		return Real.filt(Real(x)*Real(180)/Real.pi())

	def rad(x) :
		return Real.filt(Real(x)*Real.pi()/Real(180))

	def sqrt(x) :
		x = Real(x)
		return Decimal.sqrt(x)

	def __repr__(self):
		return f"Real({self})"

	def filt(self):
		(sign, ints, low) = self.as_tuple()
		length = len(ints)
		if self == 0 :
			low = 0
			ints = (0,)
		if isinstance(low,int) :
			while low<0 and ints[len(ints)-1] == 0 :
				ints = ints[:len(ints)-1]
				low += 1
			if len(ints)==1 and ints[0] == 0 :
				low = 0
		return Real( (sign, ints, low) )
		
	def is_int(self) :
		if isinstance(self,int) : return True
		x = Real.filt(self)
		(sign, ints, low) = x.as_tuple()
		if isinstance(low,int) and low >= 0 : return True
		return False

	#def __str__(self):
	#	(sign, ints, exp) = self.as_tuple()
	#	if  not isinstance(exp,int) :
	#		return Decimal.str(self)
	#	txt = Complex('-' if sign>0 else ''
	#	length = len(ints)
	#	idx = 0
	#	for num in ints :
	#		if length + exp == idx :
	#			txt += '.'
	#		txt += str(num)
	#		idx += 1
	#	txt += ")"

"""
The high precision number class for Complex numbers
"""
class Complex(object):

	def testprint(self,text) :
		print( self , text)
		#if isinstance(self,Complex) :
		#	(signreal, intsreal, expreal) = self.real.as_tuple()
		#	lenreal = len(intsreal)
		#	(signimag, intsimag, expimag) = self.imag.as_tuple()
		#	lenimag = len(intsimag)
		#	print('(',signreal, intsreal, expreal,')','(',signimag, intsimag, expimag,')',text)
		#elif isinstance(self,Real) :
		#	(sign, ints, exp) = self.as_tuple()
		#	lenreal = len(ints)
		#	print('(',sign, ints, exp,')',text)
		return

	def __init__(self,real,imag=0):
		if type(real) == type(2.0+2.j) :
			self.real = Real(real.real)
			self.imag = Real(real.imag)
		elif type(real) == type(Complex(2.0j)) :
			self.real = real.real
			self.imag = real.imag
		else :
			self.real = Real(real)
			self.imag = Real(imag)

	def pi() :
		return Complex(Real.pi())

	def e():
		return Complex(Real.e())

	def j():
		return Complex(1j)

	def conj(self): return Complex(self.real,-self.imag)

	def re(self): return self.real

	def im(self): return self.imag

	def floor(self):
		return Complex(self.real.to_integral(), self.imag.to_integral())

	def ceil(self): # cut?
		return Complex(math.ceil(self.real),math.ceil(self.imag))

	def round(self,prec=0):
		return Complex(round(self.real,prec),round(self.imag,prec))

	def exp(self):
		#self.testprint('exp')
		if self.imag == 0:
			return Complex(exp(self.real))
		cos = Real.cos(self.imag)
		sin = Real.sin(self.imag)
		with add_prec():
			e = Complex( cos, sin ) * self.real.exp()
		return +e

	def phase(self):
		#self.testprint('phase')
		(abso,phase) = self.polar()
		return phase

	def polar(self): #Converts to polar
		abso = abs(self)
		sin = self.imag/abso
		phase = Real.asin(sin)
		if self.real < 0 :
			if phase < 0 :
				phase = Real.pi()/2 - phase
			else :
				phase = phase - Real.pi()/2
		#Complex.testprint(Real.grad(phase),'polar arg')
		return (abso, phase)

	def rect(self,r,phi):
		return Complex(Real.cos(phi)*r , Real.sin(phi)*r)

	def cos(self):
		if self.imag == 0:
			return Complex(Real.cos(self.real))
		a = j() * self
		return (Complex.exp(-a) + Complex.exp(a))*0.5

	def sin(self):
		if self.imag == 0:
			return Complex(Real.sin(self.real))
		a = j() * self
		return j() * (Complex.exp(-a) - Complex.exp(a))*0.5

	def tan(self):
		try: return self.sin() / self.cos()
		except: return Complex("Inf")

	def acos(self):
		if self.imag == 0:
			return Complex(math.aacos(self))
		A = (((1 + self.real)**2 + self.imag**2).sqrt() - ((1 - self.real).sqrt() + self.imag**2).sqrt()) / 2
		B = (((1 + self.real)**2 + self.imag**2).sqrt() + ((1 - self.real).sqrt() + self.imag**2).sqrt()) / 2
		return Complex(A.acos(), -(B+(B*B - 1).sqrt()).ln() )

	def asin(self):
		if self.imag == 0:
			return Complex(math.aasin(self))
		return - ( ( 1 - self*self).sqrt() + j() * self ).ln() * j()

	def atan(self):
		#self.testprint('atan')
		result = (Complex(0,1)+self) / (Complex(0,1) - self)
		return Complex(0,0.5) * result.ln()

	def sinh(self):
		result = (self.exp() - (-self).exp()) / 2
		result.imag = -result.imag
		return result

	def cosh(self):
		result = (self.exp() + (-self).exp()) / 2
		result.imag = -result.imag
		return result

	def tanh(self):
		return self.sinh()/self.cosh()

	def acosh(self):
		result = self + (self*self - 1)**0.5
		return result.ln()

	def asinh(self):
		result = self + (self*self + 1)**0.5
		return result.ln()

	def atanh(self):
		result = (1+self).ln() - (1-self).ln()
		return result * 0.5

	def ln(self):
		getcontext().prec += 2
		(r,fi) = Complex.polar(self)
		l = Complex(Real.ln(r), fi)
		getcontext().prec -= 2
		return +l

	def __pow__(self,other):
		self = +self
		#self.testprint('**base')
		#Complex.testprint(other,'**exponent')
		if isinstance(other,(int,float,Real,Decimal)) or other.imag==0 :
			# real or int exponent
			exp = +Real(other)
			if self.imag==0 :
				# real base, real or int exponent
				return Complex(self.real**exp,0)
			# complex base, real or int exponent
			a = exp.as_tuple()
			if a[2]>=0 and exp >= 0 : # positive int exponent
				expi = int(exp) # limited exponent !!!
				if expi == 0 :
					return Complex(1)
				if expi == 1 :
					return self
				result = Complex(0)
				n = 0
				while n <= expi :
					tag = binom( expi, n )
					i = 0
					while i < n :
						tag *= self.imag 
						i += 1
					while i < expi :
						tag *= self.real
						i += 1
					match n%4:
						case 0:
							result += tag
						case 1:
							result += tag * j()
						case 2:
							result -= tag
						case 3:
							result -= tag * j()
					n += 1
				return result
			# complex base, real or negative int exponent
			exp = Complex(other) if isinstance(other,complex) else other
		else : #Complex exponent
			exp = Complex(other) if isinstance(other,complex) else other
		p = Complex.exp( Complex.ln(self) * exp)
		return +p

	def __rpow__(self,other):
		return Complex(other) ** self

	#def log(self,x=None):
	#	if not x: x = Real.e()
	#	if self.imag == 0:
	#		try: return Complex(math.log(self.real,x))
	#		except: pass
	#	p = self.polar()
	#	return Complex(math.log(p[0],x), p[1])

	#def log10(self):
	#	if self.imag == 0:
	#		return Complex(math.log10(self.real))
	#	return self.ln() / Complex(10).ln()

	def __str__(self):
		if self.imag == 0 and self.real == 0 : return '0'
		r = str(int(self.real)) if self.real.exp==0 else str(self.real)
		i = str(int(self.imag)) if self.imag.exp==0 else str(self.imag)
		if self.imag == 0: return r
		if self.real == 0: return i+'j'
		c = '+' if self.imag.as_tuple()[0]==0 else ''
		return r+c+i+'j'

	def __repr__(self):
		"""Represents the number as an instance of Decimal."""
		return "Complex('%s')" % str(self)

	def __abs__(self):
		return (self.real**2 + self.imag**2).sqrt()

	def abs(self):
		return (self.real**2 + self.imag**2).sqrt()

	def __add__(self,other):
		if type(other) == type(self) :
			return Complex(self.real+other.real, self.imag+other.imag)
		else:
			return self + Complex(other)

	def __radd__(self,other):
		return Complex(other) + self

	def __sub__(self,other):
		if type(other) == type(self) :
			return Complex(self.real-other.real, self.imag-other.imag)
		else:
			return self - Complex(other)

	def __rsub__(self,other):
		return Complex(other) - self

	def __mul__(self,other):
		if type(other) == type(self) :
			if self.imag == 0 and other.imag == 0:
				return Complex(self.real*other.real)
			return Complex( self.real*other.real - self.imag*other.imag, self.real*other.imag + self.imag*other.real )
		else:
			return self * Complex(other)

	def __rmul__(self,other):
		return Complex(other) * self

	def __floordiv__(self,other):
		return Complex(self.real // other.real)

	def __rfloordiv__(self,other):
		return Complex(other) // self

	def __truediv__(self,other):
		if type(other) == type(self) :
			if self.imag == 0 and other.imag == 0:
				return Complex(self.real/other.real)
			with add_prec():
				a = self.real*other.real + self.imag*other.imag
				b = self.imag*other.real - self.real*other.imag
				c = 1 / (other.real*other.real + other.imag*other.imag)
				d = Complex(a*c,b*c)
			return +d
		else:
			return self / Complex(other)

	def __rtruediv__(self,other):
		return Complex(other) / self

	def __neg__(self):
		return Complex(-self.real,-self.imag)

	def __pos__(self):
		''' normalize numbers together '''
		(signreal, intsreal, lowreal) = self.real.as_tuple()
		(signimag, intsimag, lowimag) = self.imag.as_tuple()
		if self.real == 0 and self.imag == 0:
			if len(intsreal)==1 and intsreal[0] == 0 :
				lowreal = 0
				self.real = Real( (signreal, intsreal, lowreal) )
			if len(intsimag)==1 and intsimag[0] == 0 :
				lowimag = 0
				self.imag = Real( (signimag, intsimag, lowimag) )
			return self
		if self.real == 0 or type(self.real.as_tuple()[2])==type('a'):
			if len(intsreal)==1 and intsreal[0] == 0 :
				lowreal = 0
				self.real = Real( (signreal, intsreal, lowreal) )
			self.imag = +self.imag;
			(signimag, intsimag, lowimag) = self.imag.as_tuple()
			mod = False
			while lowimag<0 and intsimag[-1] == 0 :
				intsimag = intsimag[:-1]
				lowimag += 1
				mod = True
			if len(intsimag)==1 and intsimag[0] == 0 :
				lowimag = 0
				mod = True
			if mod :
				self.imag = Real( (signimag, intsimag, lowimag) )
			return self
		if self.imag == 0 or type(self.imag.as_tuple()[2])==type('a'):
			self.real = +self.real;
			(signreal, intsreal, lowreal) = self.real.as_tuple()
			mod = False
			while lowreal<0 and intsreal[-1] == 0 :
				intsreal = intsreal[:-1]
				lowreal += 1
				mod = True
			if len(intsreal)==1 and intsreal[0] == 0 :
				lowreal = 0
				mod = True
			if mod :
				self.real = Real( (signreal, intsreal, lowreal) )
			return self

		# tuple (sign, ints_tuple, exponent) -- (-1)**sign * ints * 10**exponent
		heightreal = lowreal + len(intsreal)
		heightimag = lowimag + len(intsimag)

		heightcommon = max(heightreal, heightimag)
		lowcommon = heightcommon - getcontext().prec+1
		
		if heightreal > lowcommon :
			mod = False
			if heightreal-len(intsreal) < lowcommon :
				intsreal = intsreal[:heightreal-lowcommon]
				lowreal = heightreal - len(intsreal)
				mod = True
			while lowreal<0 and intsreal[-1] == 0 :
				intsreal = intsreal[:-1]
				lowreal += 1
				mod = True
			if len(intsreal)==1 and intsreal[0] == 0 :
				lowreal = 0
				mod = True
			if mod :
				self.real = Real( (signreal, intsreal, lowreal) )
		else :
			self.real = Real(0)

		if heightimag > lowcommon :
			mod = False
			if heightimag-len(intsimag) < lowcommon :
				intsimag = intsimag[:heightimag-lowcommon]
				lowimag = heightimag - len(intsimag)
				mod = True
			while lowimag<0 and intsimag[-1] == 0 :
				intsimag = intsimag[:-1]
				lowimag += 1
				mod = True
			if len(intsimag)==1 and intsimag[0] == 0 :
				lowimag = 0
				mod = True
			if mod :
				self.imag = Real( (signimag, intsimag, lowimag) )
		else :
			self.imag = Real(0)

		return self

	def __inverse__(self):
		return 1/self

	def __mod__(self,other):#a%b = a + b * ciel(-a/b) 
		return self + other * ((-self/other).ceil())

	def __complex__(self): 
		return complex(float(self.real),float(self.imag))
	def __int__(self):
		return int(self.real)
	def __float__(self):
		return float(self.real)

	def __lt__(self,other): #>
		if self.imag != 0 or other.imag != 0:
			raise errors.ComparasionError("Complex comparasion is not supported")
		return self.real < other.real

	def __le__(self,other): #>=
		if self.imag != 0 or other.imag != 0:
			raise errors.ComparasionError("Complex comparasion is not supported")
		return self.real <= other.real

	def __eq__(self,other): #==
		if self.real == other.real and self.imag == other.imag:
			return True
		return False

	def __ne__(self,other): #!=
		return not self.__eq__(other)

	def __gt__(self,other): #<
		if self.imag != 0 or other.imag != 0:
			raise errors.ComparasionError("Complex comparasion is not supported")
		return self.real > other.real

	def __ge__(self,other): #<=
		if self.imag != 0 or other.imag != 0:
			raise errors.ComparasionError("Complex comparasion is not supported")
		return self.real >= other.real

	#Some things reimplemented from decimal class for complex CDecimals
	#================================================================
	def copy_abs(self): return abs(self)

	def copy_negate(self): return -self

	def copy_sign(self,other):
		if other < 0:
			return -abs(self)
		return abs(self)

	def is_finite(self):
		return self.real.is_finite() and self.imag.is_finite()

	def is_infinite(self):
		return not self.is_finite()

	def is_nan(self):
		return self.real.is_nan() or self.imag.is_nan()

	def is_signed(self):
		return self < Complex(0)

	def is_zero(self):
		return self.real.is_zero() and self.imag.is_zero()

	def radix(self): return Decimal(10)
	def sqrt(self): 
		if self.imag.is_zero():
			return Complex(self.real.sqrt()) if self.real>0 else Complex(0,abs(self).real.sqrt())
		return self ** 0.5

	def to_eng_string(self): 
		returned = self.real.to_eng_string() + " + " + self.imag.to_eng_string + "j"
		return returned.replace(" + -"," - ")


class ComparasionError(Exception):
    pass




if __name__ == '__main__' : # short module test

	def ex(txt) :
		txt1 = txt.replace("'",'"')
		txt = "print('" + txt1 + " = '.replace('()',''), " + txt + ")"
		exec(txt)

	ex(" Real(1).is_int() ")
	ex(" Real(1000).is_int() ")
	ex(" Real('1e500').is_int() ")

	getcontext().prec = 30
	pi = Real.pi()
	sq = Real.sqrt
	sq32 = sq(3)/2
	print(getcontext().prec, ' digits')
	ex(" Real.atan2( 0,    1   ).grad()  ")
	ex(" Real.atan2( 0.5,  sq32).grad()  ")
	ex(" Real.atan2( 1,    1   ).grad()  ")
	ex(" Real.atan2( sq32, 0.5 ).grad()  ")
	ex(" Real.atan2( 1,    0   ).grad()  ")
	ex(" Real.atan2( sq32,-0.5 ).grad()  ")
	ex(" Real.atan2( 1,   -1   ).grad()  ")
	ex(" Real.atan2( 0.5, -sq32).grad()  ")
	ex(" Real.atan2( 0,   -1   ).grad()  ")
	ex(" Real.atan2(-0.5, -sq32).grad()  ")
	ex(" Real.atan2(-1,   -1   ).grad()  ")
	ex(" Real.atan2(-sq32,-0.5 ).grad()  ")
	ex(" Real.atan2(-1,    0   ).grad()  ")
	ex(" Real.atan2(-sq32, 0.5 ).grad()  ")
	ex(" Real.atan2(-1,    1   ).grad()  ")
	ex(" Real.atan2(-0.5,  sq32).grad()  ")

	ex(" e() ** ( j() * pi / 2) ")
	ex(" pi ")
	ex(" e()  ")
	ex(" Complex(0.7,0.75) ** 196 ")
	ex(" Complex(0.7,0.75) ** 197 ")
	ex(" Complex(0.7,0.75) ** 198 ")
	ex(" Complex(0.7,0.75) ** 199 ")
	ex(" Complex(0.7,0.75) ** 200 ")
	ex(" Complex(0.7,0.75) ** 201 ")
	ex(" Complex(0.7,0.75) ** 202 ")
	ex(" Complex(0.7,0.75) ** 203 ")
	ex(" Complex(0.7,0.75) ** 204 ")
	ex(" Complex(0.7,0.75) ** 205 ")
	ex(" Complex(0.7,0.75) ** 206 ")
	ex(" Complex(0.7,0.75) ** 207 ")
	ex(" Complex(0.7,0.75) ** 208 ")
	ex(" Complex(0.7,0.75) ** 209 ")
	ex(" Complex(0.7,0.75) ** 210 ")
	ex(" Complex(0.7,0.75) ** 211 ")
	ex(" Complex(0.7,0.75) ** 212 ")

	getcontext().prec = 35
	print(getcontext().prec, ' digits')
	pi = Real.pi()
	sq = Real.sqrt
	ex(" e() ** ( j() * pi ) ")
	ex(" (10 + 5*j() ) ** (2 + j() * pi ) ")

	getcontext().prec = 45
	print(getcontext().prec, ' digits')
	ex(" e() ** ( j() * pi * 3 / 2) ")
	ex(" e() ** ( j() * pi * 2) ")
	ex(" Complex.exp(j() * pi )  ")
	ex(" Complex.asin(Complex(Real.pi(),1)) ")
	ex(" Complex(1,1).asin( ) ")
	ex(" Complex(50,0).sqrt( ) ")
	ex(" Complex(-50,0).sqrt( ) ")
	ex(" Complex(-50,50).sqrt( ) ")
	ex(" binom(4,2)")
	ex(" binom(8,3)")
	ex(" binom(10,5)")
	ex(" binom(20,15)")
	ex(" binom(32+0j,15)")
	ex(" binom(32-3j,15)")
	ex(" binom(36+0j,18)")
	ex(" binom(56+0j,27)")
	ex(" binom(64,32)")
	ex(" binom(128,64)")
	ex(" binom(300+0j,150)")
	ex(" binom(599,299)")
	ex(" binom(Real(600),300)")
	ex(" Complex(30-40j).abs()")
	ex(" Complex(12-5j).abs()")
	print()

	getcontext().prec = 40
	print(getcontext().prec, ' digits')
	pi = Real.pi()
	sq = Real.sqrt
	angles  = ( 0,           pi/12,    pi/6,    pi/4,    pi/3,         5*pi/12,        pi/2)
	sinus   = ( 0, (sq(6)-sq(2))/4,     1/2, sq(2)/2, sq(3)/2, (sq(6)+sq(2))/4,           1)
	tangens = ( 0,         2-sq(3), sq(3)/3,       1,   sq(3),         2+sq(3), Real('inf'))
	i = 0
	while i < len(angles) :
		angle = Real(angles[i])

		sin_elm = Real(sinus[i])
		print(sin_elm,' good sinus at:',angle,'(',angle.grad(),'°)')
		sin_calc = Real.sin(angle)
		sinus_err = Real.filt(sin_calc - sin_elm)
		print(sin_calc,' calculated sin,  error:', sinus_err)

		tan_elm = Real(tangens[i])
		print(tan_elm,' good tangens')
		tan_calc = Real.tan(angle)
		tangens_err = Real.filt(tan_calc - tan_elm)
		print(tan_calc,' calculated tan,  error:', tangens_err)

		asin_calc = Real.asin(sin_elm)
		asinus_err = Real.filt(asin_calc - angle)
		print(asin_calc,' calculated asin,  error:', asinus_err)

		atan_calc = Real.atan(tan_elm)
		atangens_err = Real.filt(atan_calc - angle)
		print(atan_calc,' calculated atan,  error:', atangens_err)

		print()
		i += 1

	ex(" pi ")
	ex(" e()  ")
