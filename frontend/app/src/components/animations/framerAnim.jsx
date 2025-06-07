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
};;

