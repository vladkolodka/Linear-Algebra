<?php
namespace vladkolodka\linearAlgebra;

class SingularValueDecomposition {
    /**
     *    Internal storage of U.
     * @var array
     */
    private $U = array();

    /**
     * Internal storage of V.
     * @var array
     */
    private $V = array();

    private $W = array();
    /**
     * Internal storage of singular values.
     * @var array
     */
    private $s = array();

    /**
     * Row dimension.
     * @var int
     */
    private $m;

    /**
     * Column dimension.
     * @var int
     */
    private $n;

    private $K;
    private $Rank;

    /**
     * @param $Matrix Matrix
     */
    public function __construct($Matrix) {

        $this->m = $Matrix->getRowDimension();
        $this->n = $Matrix->getColumnDimension();

        $this->U = $Matrix->matrix;
        $this->V = $Matrix->getMatrix($this->n, $this->n)->matrix;

        $eps = 2.22045e-016;

        // Decompose Phase

        // Householder reduction to bi-diagonal form.
        $g = $scale = $aNorm = 0.0;
        $l = 0;
        $rv1 = [];
        for ($i = 0; $i < $this->n; $i++) {
            $l = $i + 2;
            $rv1[$i] = $scale * $g;
            $g = $s = $scale = 0.0;
            if ($i < $this->m) {
                for ($k = $i; $k < $this->m; $k++) $scale += abs($this->U[$k][$i]);
                if ($scale != 0.0) {
                    for ($k = $i; $k < $this->m; $k++) {
                        $this->U[$k][$i] /= $scale;
                        $s += $this->U[$k][$i] * $this->U[$k][$i];
                    }
                    $f = $this->U[$i][$i];
                    $g = -Matrix::sameSign(sqrt($s), $f);
                    $h = $f * $g - $s;
                    $this->U[$i][$i] = $f - $g;
                    for ($j = $l - 1; $j < $this->n; $j++) {
                        for ($s = 0.0, $k = $i; $k < $this->m; $k++) $s += $this->U[$k][$i] * $this->U[$k][$j];
                        $f = $s / $h;
                        for ($k = $i; $k < $this->m; $k++) $this->U[$k][$j] += $f * $this->U[$k][$i];
                    }
                    for ($k = $i; $k < $this->m; $k++) $this->U[$k][$i] *= $scale;
                }
            }
            $this->W[$i] = $scale * $g;
            $g = $s = $scale = 0.0;
            if ($i + 1 <= $this->m && $i + 1 != $this->n) {
                for ($k = $l - 1; $k < $this->n; $k++) $scale += abs($this->U[$i][$k]);
                if ($scale != 0.0) {
                    for ($k = $l - 1; $k < $this->n; $k++) {
                        $this->U[$i][$k] /= $scale;
                        $s += $this->U[$i][$k] * $this->U[$i][$k];
                    }
                    $f = $this->U[$i][$l - 1];
                    $g = -Matrix::sameSign(sqrt($s), $f);
                    $h = $f * $g - $s;
                    $this->U[$i][$l - 1] = $f - $g;
                    for ($k = $l - 1; $k < $this->n; $k++) $rv1[$k] = $this->U[$i][$k] / $h;
                    for ($j = $l - 1; $j < $this->m; $j++) {
                        for ($s = 0.0, $k = $l - 1; $k < $this->n; $k++) $s += $this->U[$j][$k] * $this->U[$i][$k];
                        for ($k = $l - 1; $k < $this->n; $k++) $this->U[$j][$k] += $s * $rv1[$k];
                    }
                    for ($k = $l - 1; $k < $this->n; $k++) $this->U[$i][$k] *= $scale;
                }
            }
            $aNorm = max($aNorm, (abs($this->W[$i]) + abs($rv1[$i])));
        }

        // Accumulation of right-hand transformations.
        for ($i = $this->n - 1; $i >= 0; $i--) {
            if ($i < $this->n - 1) {
                if ($g != 0.0) {
                    for ($j = $l; $j < $this->n; $j++) // Double division to avoid possible underflow.
                        $this->V[$j][$i] = ($this->U[$i][$j] / $this->U[$i][$l]) / $g;
                    for ($j = $l; $j < $this->n; $j++) {
                        for ($s = 0.0, $k = $l; $k < $this->n; $k++) $s += ($this->U[$i][$k] * $this->V[$k][$j]);
                        for ($k = $l; $k < $this->n; $k++) $this->V[$k][$j] += $s * $this->V[$k][$i];
                    }
                }
                for ($j = $l; $j < $this->n; $j++) $this->V[$i][$j] = $this->V[$j][$i] = 0.0;
            }
            $this->V[$i][$i] = 1.0;
            $g = $rv1[$i];
            $l = $i;
        }

        // Accumulation of left-hand transformations.
        for ($i = min($this->m, $this->n) - 1; $i >= 0; $i--) {
            $l = $i + 1;
            $g = $this->W[$i];
            for ($j = $l; $j < $this->n; $j++) $this->U[$i][$j] = 0.0;
            if ($g != 0.0) {
                $g = 1.0 / $g;
                for ($j = $l; $j < $this->n; $j++) {
                    for ($s = 0.0, $k = $l; $k < $this->m; $k++) $s += $this->U[$k][$i] * $this->U[$k][$j];
                    $f = ($s / $this->U[$i][$i]) * $g;
                    for ($k = $i; $k < $this->m; $k++) $this->U[$k][$j] += $f * $this->U[$k][$i];
                }
                for ($j = $i; $j < $this->m; $j++) $this->U[$j][$i] *= $g;
            } else {
                for ($j = $i; $j < $this->m; $j++) $this->U[$j][$i] = 0.0;
            }
            ++$this->U[$i][$i];
        }

        // Diagonalization of the bi-diagonal form
        // Loop over singular values, and over allowed iterations.
        $nm = 0;
        for ($k = $this->n - 1; $k >= 0; $k--) {
            for ($its = 0; $its < 30; $its++) {
                $flag = true;
                for ($l = $k; $l >= 0; $l--) {
                    $nm = $l - 1;
                    if ($l == 0 || abs($rv1[$l]) <= $eps * $aNorm) {
                        $flag = false;
                        break;
                    }
                    if (abs($this->W[$nm]) <= $eps * $aNorm) break;
                }
                if ($flag) {
                    $c = 0.0;  // Cancellation of rv1[l], if l > 0.
                    $s = 1.0;
                    for ($i = $l; $i < $k + 1; $i++) {
                        $f = $s * $rv1[$i];
                        $rv1[$i] = $c * $rv1[$i];
                        if (abs($f) <= $eps * $aNorm) break;
                        $g = $this->W[$i];
                        $h = Matrix::hypo($f, $g);
                        $this->W[$i] = $h;
                        $h = 1.0 / $h;
                        $c = $g * $h;
                        $s = -$f * $h;
                        for ($j = 0; $j < $this->m; $j++) {
                            $y = $this->U[$j][$nm];
                            $z = $this->U[$j][$i];
                            $this->U[$j][$nm] = $y * $c + $z * $s;
                            $this->U[$j][$i] = $z * $c - $y * $s;
                        }
                    }
                }
                $z = $this->W[$k];
                if ($l == $k) {
                    if ($z < 0.0) {
                        $this->W[$k] = -$z; // Singular value is made nonnegative.
                        for ($j = 0; $j < $this->n; $j++) $this->V[$j][$k] = -$this->V[$j][$k];
                    }
                    break;
                }
                if ($its == 29) print("no convergence in 30 svd iterations");
                $x = $this->W[$l]; // Shift from bottom 2-by-2 minor.
                $nm = $k - 1;
                $y = $this->W[$nm];
                $g = $rv1[$nm];
                $h = $rv1[$k];
                $f = (($y - $z) * ($y + $z) + ($g - $h) * ($g + $h)) / (2.0 * $h * $y);
                $g = Matrix::hypo($f, 1.0);
                $f = (($x - $z) * ($x + $z) + $h * (($y / ($f + Matrix::sameSign($g, $f))) - $h)) / $x;
                $c = $s = 1.0;
                for ($j = $l; $j <= $nm; $j++) {
                    $i = $j + 1;
                    $g = $rv1[$i];
                    $y = $this->W[$i];
                    $h = $s * $g;
                    $g = $c * $g;
                    $z = Matrix::hypo($f, $h);
                    $rv1[$j] = $z;
                    $c = $f / $z;
                    $s = $h / $z;
                    $f = $x * $c + $g * $s;
                    $g = $g * $c - $x * $s;
                    $h = $y * $s;
                    $y *= $c;
                    for ($jj = 0; $jj < $this->n; $jj++) {
                        $x = $this->V[$jj][$j];
                        $z = $this->V[$jj][$i];
                        $this->V[$jj][$j] = $x * $c + $z * $s;
                        $this->V[$jj][$i] = $z * $c - $x * $s;
                    }
                    $z = Matrix::hypo($f, $h);
                    $this->W[$j] = $z;  // Rotation can be arbitrary if z = 0.
                    if ($z) {
                        $z = 1.0 / $z;
                        $c = $f * $z;
                        $s = $h * $z;
                    }
                    $f = $c * $g + $s * $y;
                    $x = $c * $y - $s * $g;
                    for ($jj = 0; $jj < $this->m; $jj++) {
                        $y = $this->U[$jj][$j];
                        $z = $this->U[$jj][$i];
                        $this->U[$jj][$j] = $y * $c + $z * $s;
                        $this->U[$jj][$i] = $z * $c - $y * $s;
                    }
                }
                $rv1[$l] = 0.0;
                $rv1[$k] = $f;
                $this->W[$k] = $x;
            }
        }

        // Reorder Phase
        // Sort. The method is Shell's sort.
        // (The work is negligible as compared to that already done in decompose phase.)
        $inc = 1;
        do {
            $inc *= 3;
            $inc++;
        } while ($inc <= $this->n);

        $su = [];
        $sv = [];

        do {
            $inc /= 3;
            for ($i = $inc; $i < $this->n; $i++) {
                $sw = $this->W[$i];
                for ($k = 0; $k < $this->m; $k++) $su[$k] = $this->U[$k][$i];
                for ($k = 0; $k < $this->n; $k++) $sv[$k] = $this->V[$k][$i];
                $j = $i;
                while ($this->W[$j - $inc] < $sw) {
                    $this->W[$j] = $this->W[$j - $inc];
                    for ($k = 0; $k < $this->m; $k++) $this->U[$k][$j] = $this->U[$k][$j - $inc];
                    for ($k = 0; $k < $this->n; $k++) $this->V[$k][$j] = $this->V[$k][$j - $inc];
                    $j -= $inc;
                    if ($j < $inc) break;
                }
                $this->W[$j] = $sw;
                for ($k = 0; $k < $this->m; $k++) $this->U[$k][$j] = $su[$k];
                for ($k = 0; $k < $this->n; $k++) $this->V[$k][$j] = $sv[$k];
            }
        } while ($inc > 1);

        for ($k = 0; $k < $this->n; $k++) {
            $s = 0;
            for ($i = 0; $i < $this->m; $i++) if ($this->U[$i][$k] < 0.0) $s++;
            for ($j = 0; $j < $this->n; $j++) if ($this->V[$j][$k] < 0.0) $s++;
            if ($s > ($this->m + $this->n) / 2) {
                for ($i = 0; $i < $this->m; $i++) $this->U[$i][$k] = -$this->U[$i][$k];
                for ($j = 0; $j < $this->n; $j++) $this->V[$j][$k] = -$this->V[$j][$k];
            }
        }

        // calculate the rank
        $rank = 0;
        for ($i = 0; $i < count($this->W); $i++) {
            if (round($this->W[$i], 4) > 0) {
                $rank += 1;
            }
        }

        // Low-Rank Approximation
        $q = 0.9;
        $k = 0;
        $fRobA = 0;
        $fRobAk = 0;
        for ($i = 0; $i < $rank; $i++) $fRobA += $this->W[$i];
        do {
            for ($i = 0; $i <= $k; $i++) $fRobAk += $this->W[$i];
            $clt = $fRobAk / $fRobA;
            $k++;
        } while ($clt < $q);

        // prepare S matrix as n*n daigonal matrix of singular values
        for ($i = 0; $i < $this->n; $i++) {
//                for ($j = 0; $j < $n; $j++) {
//                    $this->s[$i][$j] = 0;
            $this->s[$i]/*[$i]*/ = $this->W[$i];
//                }
        }

        $this->K = $k;
        $this->Rank = $rank;
    }

    /**
     *    Return the left singular vectors
     *
     * @access public
     * @return Matrix U
     */
    public function getU() {
        return new Matrix($this->U);
    }


    /**
     *    Return the right singular vectors
     *
     * @access public
     * @return Matrix V
     */
    public function getV() {
        return (new Matrix($this->V))->transpose();
    }


    /**
     *    Return the one-dimensional array of singular values
     *
     * @access public
     * @return array diagonal of S.
     */
    public function getSingularValues() {
        return $this->s;
    }


    /**
     *    Return the diagonal matrix of singular values
     *
     * @access public
     * @return Matrix S
     */
    public function getS() {
//        return new Matrix($this->s);
        $S = array();
        for ($i = 0; $i < $this->n; ++$i) {
            for ($j = 0; $j < $this->n; ++$j) {
                $S[$i][$j] = 0.0;
            }
            $S[$i][$i] = $this->s[$i];
        }
        return new Matrix($S);
    }


    /**
     *    Two norm
     *
     * @access public
     * @return float max(S)
     */
    public function norm2() {
        return $this->s[0];
    }


    /**
     *    Two norm condition number
     *
     * @access public
     * @return float max(S)/min(S)
     */
    public function cond() {
        return $this->s[0] / $this->s[min($this->m, $this->n) - 1];
    }


    /**
     *    Effective numerical matrix rank
     *
     * @access public
     * @return Number of non-negligible singular values.
     */
    public function rank() {
        $eps = pow(2.0, -52.0);
        $tol = max($this->m, $this->n) * $this->s[0] * $eps;
        $r = 0;
        for ($i = 0; $i < count($this->s); ++$i) {
            if ($this->s[$i] > $tol) {
                ++$r;
            }
        }
        return $r;
    }

    /**
     * @return int
     */
    public function getK() {
        return $this->K;
    }

    /**
     * @return int
     */
    public function getRank() {
        return $this->Rank;
    }
}