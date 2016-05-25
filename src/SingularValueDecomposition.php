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

    /**
     *    Construct the singular value decomposition
     *
     *    Derived from LINPACK code.
     *
     * @param $Arg Matrix Rectangular matrix
     */
    public function __construct($Arg) {
        // Initialize.
//        $A = $Arg->getArrayCopy();
        $A = $Arg->getArray();
        $this->m = $Arg->getRowDimension();
        $this->n = $Arg->getColumnDimension();
        $nu = min($this->m, $this->n);
        $e = array();
        $work = array();
        $wantU = true;
        $wantV = true;
        $nct = min($this->m - 1, $this->n);
        $nrt = max(0, min($this->n - 2, $this->m));

        // Reduce A to bi-diagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.
        for ($k = 0; $k < max($nct, $nrt); ++$k) {
            if ($k < $nct) {
                // Compute the transformation for the k-th column and
                // place the k-th diagonal in s[$k].
                // Compute 2-norm of k-th column without under/overflow.
                $this->s[$k] = 0;
                for ($i = $k; $i < $this->m; ++$i) {
                    $this->s[$k] = Matrix::hypo($this->s[$k], $A[$i][$k]);
                }
                if ($this->s[$k] != 0.0) {
                    if ($A[$k][$k] < 0.0) {
                        $this->s[$k] = -$this->s[$k];
                    }
                    for ($i = $k; $i < $this->m; ++$i) {
                        $A[$i][$k] /= $this->s[$k];
                    }
                    $A[$k][$k] += 1.0;
                }
                $this->s[$k] = -$this->s[$k];
            }

            for ($j = $k + 1; $j < $this->n; ++$j) {
                if (($k < $nct) & ($this->s[$k] != 0.0)) {
                    // Apply the transformation.
                    $t = 0;
                    for ($i = $k; $i < $this->m; ++$i) {
                        $t += $A[$i][$k] * $A[$i][$j];
                    }
                    $t = -$t / $A[$k][$k];
                    for ($i = $k; $i < $this->m; ++$i) {
                        $A[$i][$j] += $t * $A[$i][$k];
                    }
                    // Place the k-th row of A into e for the
                    // subsequent calculation of the row transformation.
                    $e[$j] = $A[$k][$j];
                }
            }

            if ($wantU and ($k < $nct)) {
                // Place the transformation in U for subsequent back
                // multiplication.
                for ($i = $k; $i < $this->m; ++$i) {
                    $this->U[$i][$k] = $A[$i][$k];
                }
            }

            if ($k < $nrt) {
                // Compute the k-th row transformation and place the
                // k-th super-diagonal in e[$k].
                // Compute 2-norm without under/overflow.
                $e[$k] = 0;
                for ($i = $k + 1; $i < $this->n; ++$i) {
                    $e[$k] = Matrix::hypo($e[$k], $e[$i]);
                }
                if ($e[$k] != 0.0) {
                    if ($e[$k + 1] < 0.0) {
                        $e[$k] = -$e[$k];
                    }
                    for ($i = $k + 1; $i < $this->n; ++$i) {
                        $e[$i] /= $e[$k];
                    }
                    $e[$k + 1] += 1.0;
                }
                $e[$k] = -$e[$k];
                if (($k + 1 < $this->m) and ($e[$k] != 0.0)) {
                    // Apply the transformation.
                    for ($i = $k + 1; $i < $this->m; ++$i) {
                        $work[$i] = 0.0;
                    }
                    for ($j = $k + 1; $j < $this->n; ++$j) {
                        for ($i = $k + 1; $i < $this->m; ++$i) {
                            $work[$i] += $e[$j] * $A[$i][$j];
                        }
                    }
                    for ($j = $k + 1; $j < $this->n; ++$j) {
                        $t = -$e[$j] / $e[$k + 1];
                        for ($i = $k + 1; $i < $this->m; ++$i) {
                            $A[$i][$j] += $t * $work[$i];
                        }
                    }
                }
                if ($wantV) {
                    // Place the transformation in V for subsequent
                    // back multiplication.
                    for ($i = $k + 1; $i < $this->n; ++$i) {
                        $this->V[$i][$k] = $e[$i];
                    }
                }
            }
        }

        // Set up the final bi-diagonal matrix or order p.
        $p = min($this->n, $this->m + 1);
        if ($nct < $this->n) {
            $this->s[$nct] = $A[$nct][$nct];
        }
        if ($this->m < $p) {
            $this->s[$p - 1] = 0.0;
        }
        if ($nrt + 1 < $p) {
            $e[$nrt] = $A[$nrt][$p - 1];
        }
        $e[$p - 1] = 0.0;
        // If required, generate U.
        if ($wantU) {
            for ($j = $nct; $j < $nu; ++$j) {
                for ($i = 0; $i < $this->m; ++$i) {
                    $this->U[$i][$j] = 0.0;
                }
                $this->U[$j][$j] = 1.0;
            }
            for ($k = $nct - 1; $k >= 0; --$k) {
                if ($this->s[$k] != 0.0) {
                    for ($j = $k + 1; $j < $nu; ++$j) {
                        $t = 0;
                        for ($i = $k; $i < $this->m; ++$i) {
                            $t += $this->U[$i][$k] * $this->U[$i][$j];
                        }
                        $t = -$t / $this->U[$k][$k];
                        for ($i = $k; $i < $this->m; ++$i) {
                            $this->U[$i][$j] += $t * $this->U[$i][$k];
                        }
                    }
                    for ($i = $k; $i < $this->m; ++$i) {
                        $this->U[$i][$k] = -$this->U[$i][$k];
                    }
                    $this->U[$k][$k] = 1.0 + $this->U[$k][$k];
                    for ($i = 0; $i < $k - 1; ++$i) {
                        $this->U[$i][$k] = 0.0;
                    }
                } else {
                    for ($i = 0; $i < $this->m; ++$i) {
                        $this->U[$i][$k] = 0.0;
                    }
                    $this->U[$k][$k] = 1.0;
                }
            }
        }

        // If required, generate V.
        if ($wantV) {
            for ($k = $this->n - 1; $k >= 0; --$k) {
                if (($k < $nrt) and ($e[$k] != 0.0)) {
                    for ($j = $k + 1; $j < $nu; ++$j) {
                        $t = 0;
                        for ($i = $k + 1; $i < $this->n; ++$i) {
                            $t += $this->V[$i][$k] * $this->V[$i][$j];
                        }
                        $t = -$t / $this->V[$k + 1][$k];
                        for ($i = $k + 1; $i < $this->n; ++$i) {
                            $this->V[$i][$j] += $t * $this->V[$i][$k];
                        }
                    }
                }
                for ($i = 0; $i < $this->n; ++$i) {
                    $this->V[$i][$k] = 0.0;
                }
                $this->V[$k][$k] = 1.0;
            }
        }

        // Main iteration loop for the singular values.
        $pp = $p - 1;
        $iterator = 0;
        $eps = pow(2.0, -52.0);

        while ($p > 0) {
            // Here is where a test for too many iterations would go.
            // This section of the program inspects for negligible
            // elements in the s and e arrays.  On completion the
            // variables case and k are set as follows:
            // case = 1  if s(p) and e[k-1] are negligible and k<p
            // case = 2  if s(k) is negligible and k<p
            // case = 3  if e[k-1] is negligible, k<p, and
            //           s(k), ..., s(p) are not negligible (qr step).
            // case = 4  if e(p-1) is negligible (convergence).
            for ($k = $p - 2; $k >= -1; --$k) {
                if ($k == -1) {
                    break;
                }
                if (abs($e[$k]) <= $eps * (abs($this->s[$k]) + abs($this->s[$k + 1]))) {
                    $e[$k] = 0.0;
                    break;
                }
            }
            if ($k == $p - 2) {
                $case = 4;
            } else {
                for ($ks = $p - 1; $ks >= $k; --$ks) {
                    if ($ks == $k) {
                        break;
                    }
                    $t = ($ks != $p ? abs($e[$ks]) : 0.) + ($ks != $k + 1 ? abs($e[$ks - 1]) : 0.);
                    if (abs($this->s[$ks]) <= $eps * $t) {
                        $this->s[$ks] = 0.0;
                        break;
                    }
                }
                if ($ks == $k) {
                    $case = 3;
                } elseif ($ks == $p - 1) {
                    $case = 1;
                } else {
                    $case = 2;
                    $k = $ks;
                }
            }
            ++$k;

            // Perform the task indicated by case.
            switch ($case) {
                // Deflate negligible s(p).
                case 1:
                    $f = $e[$p - 2];
                    $e[$p - 2] = 0.0;
                    for ($j = $p - 2; $j >= $k; --$j) {
                        $t = Matrix::hypo($this->s[$j], $f);
                        $cs = $this->s[$j] / $t;
                        $sn = $f / $t;
                        $this->s[$j] = $t;
                        if ($j != $k) {
                            $f = -$sn * $e[$j - 1];
                            $e[$j - 1] = $cs * $e[$j - 1];
                        }
                        if ($wantV) {
                            for ($i = 0; $i < $this->n; ++$i) {
                                $t = $cs * $this->V[$i][$j] + $sn * $this->V[$i][$p - 1];
                                $this->V[$i][$p - 1] = -$sn * $this->V[$i][$j] + $cs * $this->V[$i][$p - 1];
                                $this->V[$i][$j] = $t;
                            }
                        }
                    }
                    break;
                // Split at negligible s(k).
                case 2:
                    $f = $e[$k - 1];
                    $e[$k - 1] = 0.0;
                    for ($j = $k; $j < $p; ++$j) {
                        $t = Matrix::hypo($this->s[$j], $f);
                        $cs = $this->s[$j] / $t;
                        $sn = $f / $t;
                        $this->s[$j] = $t;
                        $f = -$sn * $e[$j];
                        $e[$j] = $cs * $e[$j];
                        if ($wantU) {
                            for ($i = 0; $i < $this->m; ++$i) {
                                $t = $cs * $this->U[$i][$j] + $sn * $this->U[$i][$k - 1];
                                $this->U[$i][$k - 1] = -$sn * $this->U[$i][$j] + $cs * $this->U[$i][$k - 1];
                                $this->U[$i][$j] = $t;
                            }
                        }
                    }
                    break;
                // Perform one qr step.
                case 3:
                    // Calculate the shift.
                    $scale = max(max(max(max(abs($this->s[$p - 1]), abs($this->s[$p - 2])), abs($e[$p - 2])), abs($this->s[$k])), abs($e[$k]));
                    $sp = $this->s[$p - 1] / $scale;
                    $spm1 = $this->s[$p - 2] / $scale;
                    $epm1 = $e[$p - 2] / $scale;
                    $sk = $this->s[$k] / $scale;
                    $ek = $e[$k] / $scale;
                    $b = (($spm1 + $sp) * ($spm1 - $sp) + $epm1 * $epm1) / 2.0;
                    $c = ($sp * $epm1) * ($sp * $epm1);
                    $shift = 0.0;
                    if (($b != 0.0) || ($c != 0.0)) {
                        $shift = sqrt($b * $b + $c);
                        if ($b < 0.0) {
                            $shift = -$shift;
                        }
                        $shift = $c / ($b + $shift);
                    }
                    $f = ($sk + $sp) * ($sk - $sp) + $shift;
                    $g = $sk * $ek;
                    // Chase zeros.
                    for ($j = $k; $j < $p - 1; ++$j) {
                        $t = Matrix::hypo($f, $g);
                        $cs = $f / $t;
                        $sn = $g / $t;
                        if ($j != $k) {
                            $e[$j - 1] = $t;
                        }
                        $f = $cs * $this->s[$j] + $sn * $e[$j];
                        $e[$j] = $cs * $e[$j] - $sn * $this->s[$j];
                        $g = $sn * $this->s[$j + 1];
                        $this->s[$j + 1] = $cs * $this->s[$j + 1];
                        if ($wantV) {
                            for ($i = 0; $i < $this->n; ++$i) {
                                $t = $cs * $this->V[$i][$j] + $sn * $this->V[$i][$j + 1];
                                $this->V[$i][$j + 1] = -$sn * $this->V[$i][$j] + $cs * $this->V[$i][$j + 1];
                                $this->V[$i][$j] = $t;
                            }
                        }
                        $t = Matrix::hypo($f, $g);
                        $cs = $f / $t;
                        $sn = $g / $t;
                        $this->s[$j] = $t;
                        $f = $cs * $e[$j] + $sn * $this->s[$j + 1];
                        $this->s[$j + 1] = -$sn * $e[$j] + $cs * $this->s[$j + 1];
                        $g = $sn * $e[$j + 1];
                        $e[$j + 1] = $cs * $e[$j + 1];
                        if ($wantU && ($j < $this->m - 1)) {
                            for ($i = 0; $i < $this->m; ++$i) {
                                $t = $cs * $this->U[$i][$j] + $sn * $this->U[$i][$j + 1];
                                $this->U[$i][$j + 1] = -$sn * $this->U[$i][$j] + $cs * $this->U[$i][$j + 1];
                                $this->U[$i][$j] = $t;
                            }
                        }
                    }
                    $e[$p - 2] = $f;
                    $iterator = $iterator + 1;
                    break;
                // Convergence.
                case 4:
                    // Make the singular values positive.
                    if ($this->s[$k] <= 0.0) {
                        $this->s[$k] = ($this->s[$k] < 0.0 ? -$this->s[$k] : 0.0);
                        if ($wantV) {
                            for ($i = 0; $i <= $pp; ++$i) {
                                $this->V[$i][$k] = -$this->V[$i][$k];
                            }
                        }
                    }
                    // Order the singular values.
                    while ($k < $pp) {
                        if ($this->s[$k] >= $this->s[$k + 1]) {
                            break;
                        }
                        $t = $this->s[$k];
                        $this->s[$k] = $this->s[$k + 1];
                        $this->s[$k + 1] = $t;
                        if ($wantV and ($k < $this->n - 1)) {
                            for ($i = 0; $i < $this->n; ++$i) {
                                $t = $this->V[$i][$k + 1];
                                $this->V[$i][$k + 1] = $this->V[$i][$k];
                                $this->V[$i][$k] = $t;
                            }
                        }
                        if ($wantU and ($k < $this->m - 1)) {
                            for ($i = 0; $i < $this->m; ++$i) {
                                $t = $this->U[$i][$k + 1];
                                $this->U[$i][$k + 1] = $this->U[$i][$k];
                                $this->U[$i][$k] = $t;
                            }
                        }
                        ++$k;
                    }
                    $iterator = 0;
                    --$p;
                    break;
            }
        }

    }

    /**
     *    Return the left singular vectors
     *
     * @access public
     * @return Matrix U
     */
    public function getU() {
        return new Matrix($this->U, $this->m);
    }


    /**
     *    Return the right singular vectors
     *
     * @access public
     * @return Matrix V
     */
    public function getV() {
        return new Matrix($this->V);
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
}