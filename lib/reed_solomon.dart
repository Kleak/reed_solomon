library reed_solomon;

import "package:galois_field/galois_field.dart";

List<int> rsCorrectMessage(List<int> message_in, int nsym) {
  List<int> message_out = new List.from(message_in);
  List<int> erase_pos = [];
  for (int i = 0; i < message_out.length; i++) {
    if (message_out[i] < 0) {
      message_out[i] = 0;
      erase_pos.add(i);
    }
  }
  if (erase_pos.length > nsym) return null;
  List<int> synd = _rsCalculateSyndrome(message_out, nsym);
  if (_max(synd) == 0) return message_out;
  List<int> fsynd = _rsForneySyndrome(synd, erase_pos, message_out.length);
  List<int> err_polynomial = _rsGeneratorErrorPolynomial(fsynd);
  List<int> err_pos = _rsFindErrors(err_polynomial, message_out.length);
  if (err_pos == null) return null;
  message_out = _rsCorrectErrata(message_out, synd, erase_pos..addAll(err_pos));
  synd = _rsCalculateSyndrome(message_out, nsym);
  if (_max(synd) > 0) return null;
  return message_out;
}

/**
 * Reed-Solomon main encoding function, using polynomial division (algorithm Extended Synthetic Division)
 */
List<int> rsEncodeMessage(List<int> message_in, int nsym) {
  List<int> gen = _rsGeneratorPolynomial(nsym);
  List<int> message_out =
      new List.filled(message_in.length + gen.length - 1, 0);
  message_out.setAll(0, message_in);
  for (int i = 0; i < message_in.length; i++) {
    int coef = message_out[i];
    if (coef != 0) {
      for (int j = 1; j < gen.length; j++) {
        message_out[i + j] ^= gfMultiply(gen[j], coef);
      }
    }
  }
  message_out.setAll(0, message_in);
  return []..addAll(message_out);
}

int _max(List<int> list) {
  int r = null;
  for (int i = 1; i < list.length; i++) {
    if (list[i - 1].compareTo(r == null ? list[i] : r) >= 0) {
      r = list[i];
    }
  }
  return r;
}

/**
 * Calculate the syndromes
 */
List<int> _rsCalculateSyndrome(List<int> msg, int nsym) {
  List<int> synd = new List.filled(nsym, 0);
  for (int i = 0; i < nsym; i++) {
    synd[i] = gfPolynomialEval(msg, GF_EXP[i]);
  }
  return synd;
}

/**
 * Forney algorithm, computes the values (error magnitude) to correct the input message
 */
List<int> _rsCorrectErrata(List<int> message, List<int> synd, List<int> pos) {
  List<int> coef_pos = <int>[];
  pos.forEach((int value) => coef_pos.add(message.length - 1 - value));
  List<int> loc = _rsFindErrataLocator(coef_pos);
  List<int> reversed = new List.from(synd.sublist(0, pos.length).reversed);
  List<int> eval = _rsFindErrorEvaluator(reversed, loc, pos.length - 1);
  List<int> locprime = <int>[];
  bool skipNext = false;
  locprime.addAll(loc.skip(loc.length & 1).where((int value) {
    skipNext = !skipNext;
    return skipNext;
  }));
  pos.forEach((int value) {
    int x = GF_EXP[value + GF_LOG_SIZE - message.length];
    int y = gfPolynomialEval(eval, x);
    int z = gfPolynomialEval(locprime, gfMultiply(x, x));
    int magnitude = gfDivide(y, gfMultiply(x, z));
    message[value] ^= magnitude;
  });
  return message;
}

/**
 * Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions
 * (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx"
 * with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed
 * since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the
 * erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.
 */
List<int> _rsFindErrataLocator(List<int> e_pos, {int x: null}) {
  List<int> e_loc = [1];
  for (x in e_pos) {
    e_loc = gfPolynomialMultiply(e_loc, gfPolynomialAdd([1], [GF_EXP[x], 0]));
  }
  return e_loc;
}

/**
 * Compute the error (or erasures if you supply sigma=erasures locator polynomial, or errata) evaluator polynomial Omega
 * from the syndrome and the error/erasures/errata locator Sigma.
 */
List<int> _rsFindErrorEvaluator(List<int> synd, List<int> err_loc, int nsym) {
  List<int> remainder = gfPolynomialDivide(gfPolynomialMultiply(synd, err_loc),
      [1]..addAll(new List.filled(nsym + 1, 0)));
  return remainder;
}

/**
 * Find the roots (ie, where evaluation = zero) of error polynomial by brute-force trial, this is a sort of Chien's search
 * (but less efficient, Chien's search is a way to evaluate the polynomial such that each evaluation only takes constant time)
 */
List<int> _rsFindErrors(List<int> err_loc, int nmess) {
  int errs = err_loc.length - 1;
  List<int> err_pos = <int>[];
  for (int i = 0; i < nmess; i++) {
    if (gfPolynomialEval(err_loc, GF_EXP[(GF_LOG_SIZE - 1) - i]) == 0) {
      err_pos.add(nmess - 1 - i);
    }
  }
  if (err_pos.length != errs) {
    return null;
  }
  return err_pos;
}

/**
 * Calculating the Forney syndromes
 */
List<int> _rsForneySyndrome(List<int> synd, List<int> pos, int nmess) {
  List<int> fsynd = new List.from(synd);
  pos.forEach((int value) {
    int x = GF_EXP[nmess - 1 - value];
    for (int j = 0; j < fsynd.length - 1; j++) {
      fsynd[j] = gfMultiply(fsynd[j], x) ^ fsynd[j + 1];
    }
    fsynd.removeLast();
  });
  return fsynd;
}

/**
 * Find error locator polynomial with Berlekamp-Massey algorithm
 */
List<int> _rsGeneratorErrorPolynomial(List<int> synd) {
  List<int> err_loc = [1];
  List<int> old_loc = [1];

  for (int i = 0; i < synd.length; i++) {
    old_loc.add(0);
    int delta = synd[i];
    for (int j = 1; j < err_loc.length; j++) {
      delta ^= gfMultiply(err_loc[err_loc.length - 1 - j], synd[i - j]);
    }
    if (delta != 0) {
      if (old_loc.length > err_loc.length) {
        List<int> new_loc = gfPolynomialScale(old_loc, delta);
        old_loc = gfPolynomialScale(err_loc, gfInverse(delta));
        err_loc = new_loc;
      }
      err_loc = gfPolynomialAdd(err_loc, gfPolynomialScale(old_loc, delta));
    }
  }
  err_loc.removeWhere((int value) => value == 0);
  int errs = err_loc.length - 1;
  if (errs * 2 > synd.length) {
    return null;
  }
  return err_loc;
}

/**
 * Computes the generator polynomial for a given number of error correction symbols
 */
List<int> _rsGeneratorPolynomial(int nsym) {
  List<int> g = [1];
  for (int i = 0; i < nsym; i++) {
    g = gfPolynomialMultiply(g, [1, GF_EXP[i]]);
  }
  return g;
}
