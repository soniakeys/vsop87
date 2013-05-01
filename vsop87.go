package vsop87

import (
	"errors"
	"fmt"
	"io/ioutil"
	"math"
	"strconv"
	"strings"
)

// Body constants
const (
	Sun = iota
	Mercury
	Venus
	Earth
	Mars
	Jupiter
	Saturn
	Uranus
	Neptune
	EarthMoon // barycenter
	nBodies
)

// parallel arrays, indexed by body constants.
var (
	ext = [nBodies]string{"",
		"mer", "ven", "ear", "mar", "jup", "sat", "ura", "nep", "emb"}

	b7 = [nBodies]string{
		"",
		"MERCURY",
		"VENUS  ",
		"EARTH  ",
		"MARS   ",
		"JUPITER",
		"SATURN ",
		"URANUS ",
		"NEPTUNE",
		"EMB    ",
	}

	a0 = [nBodies]float64{
		0,
		.3871,
		.7233,
		1,
		1.5237,
		5.2026,
		9.5547,
		19.2181,
		30.1096,
		1,
	}
)

// return type
type Elliptic struct {
	A float64 // semi-major axis (au)
	L float64 // mean longitude (rd)
	K float64 // k = e*cos(pi) (rd)
	H float64 // h = e*sin(pi) (rd)
	Q float64 // q = sin(i/2)*cos(omega) (rd)
	P float64 // p = sin(i/2)*sin(omega) (rd)
	// e:     eccentricity
	// pi:    perihelion longitude
	// i:     inclination
	// omega: ascending node longitude
}

// return type
type Rectangular struct {
	Px float64 // position x (au)
	Py float64 // position y (au)
	Pz float64 // position z (au)
	Vx float64 // velocity x (au/day)
	Vy float64 // velocity y (au/day)
	Vz float64 // velocity z (au/day)
}

// return type
type Spherical struct {
	Lon  float64 // longitude (rd)
	Lat  float64 // latitude (rd)
	R    float64 // radius (au)
	VLon float64 // longitude velocity (rd/day)
	VLat float64 // latitude velocity (rd/day)
	VR   float64 // radius velocity (au/day)
}

type abc struct {
	a, b, c float64
}

type EllipticModel struct {
	t [6]float64
	c [nBodies][6][6][]abc
}

const (
	t2000 = 2451545 // J2000
	a1000 = 365250  // days per Julian millenium
)

// NewEllipticModel reads VSOP87 files and returns an object that can compute
// positions.
//
// Tdj does not have to be exact.  It is used only for determining the
// subset of coefficients needed for the requested precition.
func NewEllipticModel(path string, prec, tdj float64) (*EllipticModel, error) {
	if prec < 0 || prec > .01 {
		return nil, errors.New("Invalid precision.")
	}
	q := -math.Log10(prec + 1e-50)
	if q < 3 {
		q = 3
	}
	em := &EllipticModel{}
	var at [6]float64 // abs of powers of t
	at[0] = 1
	t := math.Abs(tdj-t2000) / a1000
	for i := 1; i < 6; i++ {
		at[i] = t * at[i-1]
	}
	var c [2047]abc
	for _, ibody := range []int{Mercury, Venus, EarthMoon, Mars,
		Jupiter, Saturn, Uranus, Neptune} {
		data, err := ioutil.ReadFile(path + "/VSOP87." + ext[ibody])
		if err != nil {
			return nil, err
		}
		lines := strings.Split(string(data), "\n")
		n := 0
		dl := 0.
		for n < len(lines) {
			line := lines[n]
			if len(line) < 132 {
				break
			}
			if iv := line[17]; iv != '0' {
				return nil, fmt.Errorf("Line %d: expected version 0, "+
					"found %c.", n+1, iv)
			}
			if bo := line[22:29]; bo != b7[ibody] {
				return nil, fmt.Errorf("Line %d: expected body %s, "+
					"found %s.", n+1, b7[ibody], bo)
			}
			ic := line[41] - '0'
			it := line[59] - '0'
			in, err := strconv.Atoi(strings.TrimSpace(line[60:67]))
			if err != nil {
				return nil, fmt.Errorf("Line %d: %v.", n+1, err)
			}
			if in == 0 {
				continue
			}
			if in > len(lines)-n {
				return nil, errors.New("Unexpected end of file.")
			}
			d0 := at[it]
			p := prec / 10 / (q - 2) / (d0 + float64(it)*dl*1e-4 + 1e-50)
			p *= a0[ibody]
			dl = d0

			n++
			cx := 0
			for _, line := range lines[n : n+in] {
				a := &c[cx]
				a.a, err =
					strconv.ParseFloat(strings.TrimSpace(line[79:97]), 64)
				if err != nil {
					goto parseError
				}
				if math.Abs(a.a) < p {
					fmt.Println("truncated")
					break
				}
				a.b, err = strconv.ParseFloat(line[98:111], 64)
				if err != nil {
					goto parseError
				}
				a.c, err =
					strconv.ParseFloat(strings.TrimSpace(line[111:131]), 64)
				if err != nil {
					goto parseError
				}
				cx++
				continue
			parseError:
				return nil, fmt.Errorf("Line%d: %v.", n+cx+1, err)
			}
			em.c[ibody][ic-1][it] = append([]abc{}, c[:cx]...)
			n += in
		}
	}
	return em, nil
}

func (em *EllipticModel) Pos(tdj float64, ibody int, r *Elliptic) {
	em.t[0] = 1
	t := (tdj - t2000) / a1000
	for i := 1; i < 6; i++ {
		em.t[i] = t * em.t[i-1]
	}
	cb := em.c[ibody]
	r.A = em.sum(cb[0])
	r.L = pmod(em.sum(cb[1]), 2*math.Pi)
	r.K = em.sum(cb[2])
	r.H = em.sum(cb[3])
	r.Q = em.sum(cb[4])
	r.P = em.sum(cb[5])
}

func (em *EllipticModel) sum(c [6][]abc) (r float64) {
	for it, t := range c {
		for _, term := range t {
			u := term.b + term.c*em.t[1]
			cu := math.Cos(u)
			r += term.a * cu * em.t[it]
		}
	}
	return
}

// copied from github.com/soniakeys/meeus/base/math.go
func pmod(x, y float64) float64 {
    r := math.Mod(x, y)
    if r < 0 {
        r += y
    }
    return r
}
