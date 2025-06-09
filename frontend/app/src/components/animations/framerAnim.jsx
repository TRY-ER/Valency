import { motion } from "framer-motion";

export const fadeInUpVariants = {
  hidden: { opacity: 0, y: 20 },
  visible: (index = 0) => ({
    opacity: 1,
    y: 0,
    transition: {
      delay: index * 0.1, // adjust delay multiplier as needed
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, y: 20, transition: { duration: 0.3 } }
};

export const fadeInRightVariants = {
  hidden: { opacity: 0, x: 20 },
  visible: (index = 0) => ({
    opacity: 1,
    x: 0,
    transition: {
      delay: index * 0.1, // adjust delay multiplier as needed
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, x: 20, transition: { duration: 0.3 } }
};

export const fadeInLeftVariants = {
  hidden: { opacity: 0, x: -5 },
  visible: (index = 0) => ({
    opacity: 1,
    x: 0,
    transition: {
      delay: index * 0.1, // adjust delay multiplier as needed
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, x: -5, transition: { duration: 0.3 } }
};

export const fadeInDownVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: (index = 0) => ({
    opacity: 1,
    y: 0,
    transition: {
      delay: index * 0.1, // adjust delay multiplier as needed
      duration: 0.3,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, y: -20, transition: { duration: 0.3 } }
};

export const fadeInUpVariantStatic = {
  hidden: { opacity: 0, y: 15 },
  visible: () => ({
    opacity: 1,
    y: 0,
    transition: {
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, y: 15, transition: { duration: 0.3 } }
}

export const fadeInRightVariantStatic = {
  hidden: { opacity: 0, x: 20 },
  visible: () => ({
    opacity: 1,
    x: 0,
    transition: {
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, x: 20, transition: { duration: 0.3 } }
};

export const fadeInStatic = {
  hidden: { opacity: 0},
  visible: () => ({
    opacity: 1,
    transition: {
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, transition: { duration: 0.3 } }
};

// New enhanced animation variants

export const scaleInVariants = {
  hidden: { opacity: 0, scale: 0.8 },
  visible: (index = 0) => ({
    opacity: 1,
    scale: 1,
    transition: {
      delay: index * 0.1,
      duration: 0.4,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, scale: 0.8, transition: { duration: 0.3 } }
};

export const slideInVariants = {
  hidden: { opacity: 0, x: -100 },
  visible: (index = 0) => ({
    opacity: 1,
    x: 0,
    transition: {
      delay: index * 0.15,
      duration: 0.6,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, x: -100, transition: { duration: 0.4 } }
};

export const bounceInVariants = {
  hidden: { opacity: 0, scale: 0.3, y: 50 },
  visible: (index = 0) => ({
    opacity: 1,
    scale: 1,
    y: 0,
    transition: {
      delay: index * 0.1,
      duration: 0.6,
      ease: "easeOut",
      type: "spring",
      damping: 10,
      stiffness: 100
    },
  }),
  exit: { opacity: 0, scale: 0.3, y: 50, transition: { duration: 0.3 } }
};

export const rotateInVariants = {
  hidden: { opacity: 0, rotate: -180, scale: 0.5 },
  visible: (index = 0) => ({
    opacity: 1,
    rotate: 0,
    scale: 1,
    transition: {
      delay: index * 0.1,
      duration: 0.8,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, rotate: 180, scale: 0.5, transition: { duration: 0.4 } }
};

export const containerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: {
      staggerChildren: 0.1,
      delayChildren: 0.2,
    },
  },
  exit: {
    opacity: 0,
    transition: {
      staggerChildren: 0.05,
      staggerDirection: -1,
    },
  }
};

export const itemVariants = {
  hidden: { opacity: 0, y: 20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: {
      duration: 0.5,
      ease: "easeOut",
    },
  },
  exit: { opacity: 0, y: 20, transition: { duration: 0.3 } }
};

export const hoverVariants = {
  hover: {
    scale: 1.05,
    y: -5,
    transition: {
      duration: 0.2,
      ease: "easeInOut",
    },
  },
  tap: {
    scale: 0.95,
    transition: {
      duration: 0.1,
    },
  }
};

export const cardVariants = {
  hidden: { opacity: 0, y: 50, scale: 0.95 },
  visible: (index = 0) => ({
    opacity: 1,
    y: 0,
    scale: 1,
    transition: {
      delay: index * 0.1,
      duration: 0.5,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, y: 50, scale: 0.95, transition: { duration: 0.3 } },
  hover: {
    y: -8,
    scale: 1.02,
    transition: {
      duration: 0.2,
      ease: "easeInOut",
    },
  }
};

export const listItemVariants = {
  hidden: { opacity: 0, x: -20 },
  visible: (index = 0) => ({
    opacity: 1,
    x: 0,
    transition: {
      delay: index * 0.05,
      duration: 0.3,
      ease: "easeOut",
    },
  }),
  exit: { opacity: 0, x: -20, transition: { duration: 0.2 } }
};

export const expandableVariants = {
  collapsed: { opacity: 0, height: 0, scale: 0.95 },
  expanded: {
    opacity: 1,
    height: "auto",
    scale: 1,
    transition: {
      duration: 0.4,
      ease: "easeOut",
      height: {
        duration: 0.4,
      },
      opacity: {
        duration: 0.25,
        delay: 0.15,
      },
    },
  },
  exit: {
    opacity: 0,
    height: 0,
    scale: 0.95,
    transition: {
      duration: 0.3,
      ease: "easeIn",
    },
  }
};

