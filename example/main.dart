
import "dart:io";

import "package:galois_field/galois_field.dart";
import "package:reed_solomon/reed_solomon.dart";

printBlock(message_in, message) {
  bool isCode = false;
  for (int i = 0; i < message.length; i++) {
    if (i != 0) {
      stdout.write(" ");
    }
    if (i >= message_in.length && isCode == false) {
      isCode = true;
      stdout.write("\n");
    }
    stdout.write("0x${message[i].toRadixString(16)}");
  }
  print("");
}

main()  {
  Stopwatch sw = new Stopwatch();
  sw.start();
  initTables();

  List<int> message_in = [
    0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
    0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec
  ];
  int k = 20;
  List<int> msg = rsEncodeMessage(message_in, k);

  printBlock(message_in, msg);
  print("");

  msg[0] = 0;

  printBlock(message_in, msg);
  print("");

  msg = rsCorrectMessage(msg, k);

  printBlock(message_in, msg);


  sw.stop();
  print("time: ${sw.elapsed}");
}
