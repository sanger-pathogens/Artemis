package uk.ac.sanger.artemis.io;

class UIEraserThread implements Runnable {
   private boolean stop;

   /**
    *@param The prompt displayed to the user
    */
   public UIEraserThread(String prompt) {
       System.out.print(prompt);
   }

   /**
    * Begin masking...display asterisks (*)
    */
   public void run () {
      stop = true;
      while (stop) {
         System.out.print("\010 ");
	     try {
	         Thread.currentThread().sleep(1);
         } catch(InterruptedException ie) {
             ie.printStackTrace();
         }
      }
   }

   /**
    * Instruct the thread to stop masking
    */
   public void stopMasking() {
      this.stop = false;
   }
}